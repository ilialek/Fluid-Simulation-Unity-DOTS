using Unity.Collections;
using Unity.Entities;
using Unity.Jobs;
using Unity.Burst;
using Unity.Mathematics;
using Unity.Transforms;

namespace SPH
{
    [RequireMatchingQueriesForUpdate]
    [UpdateInGroup(typeof(SimulationSystemGroup))]
    [UpdateBefore(typeof(TransformSystemGroup))]
    public partial struct ParticleSystem : ISystem
    {
        static readonly int3[] neighborOffsets =
        {
                    new int3(-1, -1, -1),
                    new int3(-1, -1, 0),
                    new int3(-1, -1, 1),
                    new int3(-1, 0, -1),
                    new int3(-1, 0, 0),
                    new int3(-1, 0, 1),
                    new int3(-1, 1, -1),
                    new int3(-1, 1, 0),
                    new int3(-1, 1, 1),
                    new int3(0, -1, -1),
                    new int3(0, -1, 0),
                    new int3(0, -1, 1),
                    new int3(0, 0, -1),
                    new int3(0, 0, 0),
                    new int3(0, 0, 1),
                    new int3(0, 1, -1),
                    new int3(0, 1, 0),
                    new int3(0, 1, 1),
                    new int3(1, -1, -1),
                    new int3(1, -1, 0),
                    new int3(1, -1, 1),
                    new int3(1, 0, -1),
                    new int3(1, 0, 0),
                    new int3(1, 0, 1),
                    new int3(1, 1, -1),
                    new int3(1, 1, 0),
                    new int3(1, 1, 1)
                };

        private NativeArray<float3> velocities;            // Current velocities
        private NativeArray<float3> updatedVelocities;    // Updated velocities
        private EntityQuery particleQuery;
        private EntityQuery obstacleQuery;

        public void OnCreate(ref SystemState state)
        {
            // Define the queries
            particleQuery = SystemAPI.QueryBuilder().WithAll<Particle>().WithAllRW<LocalToWorld>().Build();
            obstacleQuery = SystemAPI.QueryBuilder().WithAll<Obstacle>().Build();

            // Ensure this system updates only when particles exist
            state.RequireForUpdate(particleQuery);
        }

        [BurstCompile]
        public void OnDestroy(ref SystemState state)
        {
            if (velocities.IsCreated)
                velocities.Dispose();
            if (updatedVelocities.IsCreated)
                updatedVelocities.Dispose();
        }

        [BurstCompile]
        public void OnUpdate(ref SystemState state)
        {
            int particleCount = particleQuery.CalculateEntityCount();

            // Allocate the velocities array only once
            if (velocities.IsCreated == false)
            {
                velocities = new NativeArray<float3>(particleCount, Allocator.Persistent);
                updatedVelocities = new NativeArray<float3>(particleCount, Allocator.Persistent); // Initialize the updated velocities buffer
            }

            if (particleCount > 0)
            {
                NativeArray<LocalToWorld> particlePositions = particleQuery.ToComponentDataArray<LocalToWorld>(Allocator.TempJob);
                NativeArray<Particle> particles = particleQuery.ToComponentDataArray<Particle>(Allocator.TempJob);
                NativeArray<Obstacle> obstacles = obstacleQuery.ToComponentDataArray<Obstacle>(Allocator.TempJob);
                NativeParallelMultiHashMap<int3, int> spatialGrid = new NativeParallelMultiHashMap<int3, int>(particleCount, Allocator.TempJob);

                Particle particleSettings = particles[0];
                float dt = SystemAPI.Time.DeltaTime * particleSettings.TimeScale;

                NativeArray<float3> predictedPositions = new NativeArray<float3>(particleCount, Allocator.TempJob);

                var computeGravityAndPredictedPositionsJob = new ComputeGravityAndPredictedPositionsJob
                {
                    ParticlePositions = particlePositions,
                    PredictedPositions = predictedPositions,
                    Velocities = velocities,
                    DeltaTime = dt,
                    ParticleSettings = particleSettings,
                };
                state.Dependency = computeGravityAndPredictedPositionsJob.ScheduleParallel(particleQuery, state.Dependency);

                var populateGridJob = new PopulateSpatialGridJob
                {
                    PredictedPositions = predictedPositions,
                    CellSize = particleSettings.CellRadius,
                    SpatialGrid = spatialGrid.AsParallelWriter(),
                };
                state.Dependency = populateGridJob.Schedule(particleCount, 64, state.Dependency);

                NativeArray<float2> densities = new NativeArray<float2>(particleCount, Allocator.TempJob);
                var densityJob = new ComputeDensityJob
                {
                    PredictedPositions = predictedPositions,
                    SpatialGrid = spatialGrid,
                    ParticleSettings = particleSettings,
                    Densities = densities,
                };
                state.Dependency = densityJob.Schedule(particleCount, 64, state.Dependency);

                var calculatePressureForceJob = new CalculatePressureForceJob
                {
                    PredictedPositions = predictedPositions,
                    SpatialGrid = spatialGrid,
                    Densities = densities,
                    Velocities = velocities,
                    DeltaTime = dt,
                    ParticleSettings = particleSettings,
                };
                state.Dependency = calculatePressureForceJob.Schedule(particleCount, 64, state.Dependency);

                var calculateViscosityJob = new CalculateViscosityJob
                {
                    PredictedPositions = predictedPositions,
                    SpatialGrid = spatialGrid,
                    Velocities = velocities,
                    UpdatedVelocities = updatedVelocities,
                    DeltaTime = dt,
                    ParticleSettings = particleSettings,
                };
                state.Dependency = calculateViscosityJob.Schedule(particleCount, 64, state.Dependency);

                var resolveCollisionsApplyForcesJob = new ResolveCollisionsApplyForcesJob
                {
                    ParticlePositions = particlePositions,
                    Obstacles = obstacles,
                    DeltaTime = dt,
                    Velocities = updatedVelocities,
                    ParticleSettings = particleSettings,
                    BoxSize = new float3(35, 85, 15),
                };
                state.Dependency = resolveCollisionsApplyForcesJob.ScheduleParallel(particleQuery, state.Dependency);

                state.Dependency.Complete();

                // Swap the buffers after the jobs are done
                NativeArray<float3> temp = velocities;
                velocities = updatedVelocities;
                updatedVelocities = temp;

                // Dispose of temporary allocations
                spatialGrid.Dispose();
                particles.Dispose();
                obstacles.Dispose();
                particlePositions.Dispose();
                densities.Dispose();
                predictedPositions.Dispose();
            }
        }

        [BurstCompile]
        public partial struct ComputeGravityAndPredictedPositionsJob : IJobEntity
        {
            [ReadOnly] public NativeArray<LocalToWorld> ParticlePositions;
            [ReadOnly] public float DeltaTime;
            [ReadOnly] public Particle ParticleSettings;

            public NativeArray<float3> PredictedPositions;
            public NativeArray<float3> Velocities;

            public void Execute([EntityIndexInQuery] int index, ref Particle particle)
            {
                // Use velocities array instead of particle.Velocity
                Velocities[index] += math.down() * ParticleSettings.Gravity * DeltaTime;
                PredictedPositions[index] = ParticlePositions[index].Position + Velocities[index] / 120.0f;
            }
        }

        [BurstCompile]
        public struct PopulateSpatialGridJob : IJobParallelFor
        {
            [ReadOnly] public NativeArray<float3> PredictedPositions;
            [ReadOnly] public float CellSize;
            public NativeParallelMultiHashMap<int3, int>.ParallelWriter SpatialGrid;

            public void Execute(int index)
            {
                float3 position = PredictedPositions[index];
                int3 cell = SpatialHash(position, CellSize);
                SpatialGrid.Add(cell, index);
            }
        }

        [BurstCompile]
        public struct ComputeDensityJob : IJobParallelFor
        {
            [ReadOnly] public NativeArray<float3> PredictedPositions;
            [ReadOnly] public NativeParallelMultiHashMap<int3, int> SpatialGrid;
            [ReadOnly] public Particle ParticleSettings;

            public NativeArray<float2> Densities;

            public void Execute(int index)
            {
                float3 position = PredictedPositions[index];
                int3 currentCell = SpatialHash(position, ParticleSettings.CellRadius);

                float sqrRadius = ParticleSettings.SmoothingLength * ParticleSettings.SmoothingLength;

                float density = 0f;
                float nearDensity = 0;

                for (int i = 0; i < neighborOffsets.Length; i++)
                {
                    int3 neighborCell = currentCell + neighborOffsets[i];

                    if (SpatialGrid.TryGetFirstValue(neighborCell, out int neighborIndex, out var iterator))
                    {
                        do
                        {
                            float3 neighbourPos = PredictedPositions[neighborIndex];
                            float3 offsetToNeighbour = neighbourPos - position;
                            float sqrDstToNeighbour = math.dot(offsetToNeighbour, offsetToNeighbour);

                            if (sqrDstToNeighbour > sqrRadius) continue;

                            float dst = math.sqrt(sqrDstToNeighbour);
                            density += SpikyKernelPow2(dst, ParticleSettings.SmoothingLength);
                            nearDensity += SpikyKernelPow3(dst, ParticleSettings.SmoothingLength);

                        } while (SpatialGrid.TryGetNextValue(out neighborIndex, ref iterator));
                    }
                }

                Densities[index] = new float2(density, nearDensity);
            }
        }

        [BurstCompile]
        public struct CalculatePressureForceJob : IJobParallelFor
        {
            [ReadOnly] public NativeArray<float3> PredictedPositions;
            [ReadOnly] public NativeArray<float2> Densities;
            [ReadOnly] public NativeParallelMultiHashMap<int3, int> SpatialGrid;
            [ReadOnly] public float DeltaTime;
            [ReadOnly] public Particle ParticleSettings;

            public NativeArray<float3> Velocities;

            public void Execute(int index)
            {
                float3 position = PredictedPositions[index];
                int3 currentCell = SpatialHash(position, ParticleSettings.CellRadius);

                float density = Densities[index].x;
                float densityNear = Densities[index].y;
                float pressure = ParticleSettings.Stiffness * (density - ParticleSettings.ReferenceDensity);
                float pressureNear = densityNear * ParticleSettings.NearPressureMultiplier;

                float3 pressureForce = float3.zero;

                float sqrRadius = ParticleSettings.SmoothingLength * ParticleSettings.SmoothingLength;

                for (int i = 0; i < neighborOffsets.Length; i++)
                {
                    int3 neighborCell = currentCell + neighborOffsets[i];

                    if (SpatialGrid.TryGetFirstValue(neighborCell, out int neighborIndex, out var iterator))
                    {
                        do
                        {
                            if (neighborIndex == index) continue;

                            float3 neighbourPos = PredictedPositions[neighborIndex];
                            float3 offsetToNeighbour = neighbourPos - position;
                            float sqrDstToNeighbour = math.dot(offsetToNeighbour, offsetToNeighbour);

                            // Skip if not within radius
                            if (sqrDstToNeighbour > sqrRadius) continue;

                            float neighborDensity = Densities[neighborIndex].x;
                            float neighborNearDensity = Densities[neighborIndex].y;

                            float neighbourPressure = ParticleSettings.Stiffness * (neighborDensity - ParticleSettings.ReferenceDensity);
                            float neighbourPressureNear = neighborNearDensity * ParticleSettings.NearPressureMultiplier;

                            float sharedPressure = (pressure + neighbourPressure) / 2;
                            float sharedNearPressure = (pressureNear + neighbourPressureNear) / 2;

                            float dst = math.sqrt(sqrDstToNeighbour);
                            float3 dir = dst > 0 ? offsetToNeighbour / dst : new float3(0, 1, 0);

                            pressureForce += dir * DerivativeSpikyPow2(dst, ParticleSettings.SmoothingLength) * sharedPressure / neighborDensity;
                            pressureForce += dir * DerivativeSpikyPow3(dst, ParticleSettings.SmoothingLength) * sharedNearPressure / neighborNearDensity;

                        } while (SpatialGrid.TryGetNextValue(out neighborIndex, ref iterator));

                    }
                }

                float3 acceleration = pressureForce / density;
                // Use the velocities array to update velocity instead of particle.Velocity
                Velocities[index] += acceleration * DeltaTime;
            }
        }

        [BurstCompile]
        public struct CalculateViscosityJob : IJobParallelFor
        {
            [ReadOnly] public NativeArray<float3> PredictedPositions;
            [ReadOnly] public NativeParallelMultiHashMap<int3, int> SpatialGrid;

            [ReadOnly] public NativeArray<float3> Velocities; // Current velocities
            public NativeArray<float3> UpdatedVelocities; // Updated velocities

            [ReadOnly] public float DeltaTime;
            [ReadOnly] public Particle ParticleSettings;

            public void Execute(int index)
            {
                float3 position = PredictedPositions[index];
                int3 currentCell = SpatialHash(position, ParticleSettings.CellRadius);

                float3 viscosityForce = float3.zero;
                float sqrRadius = ParticleSettings.SmoothingLength * ParticleSettings.SmoothingLength;

                for (int i = 0; i < neighborOffsets.Length; i++)
                {
                    int3 neighborCell = currentCell + neighborOffsets[i];

                    if (SpatialGrid.TryGetFirstValue(neighborCell, out int neighborIndex, out var iterator))
                    {
                        do
                        {
                            if (neighborIndex == index) continue;

                            float3 neighbourPos = PredictedPositions[neighborIndex];
                            float3 offsetToNeighbour = neighbourPos - position;
                            float sqrDstToNeighbour = math.dot(offsetToNeighbour, offsetToNeighbour);

                            // Skip if not within radius
                            if (sqrDstToNeighbour > sqrRadius) continue;

                            float dst = math.sqrt(sqrDstToNeighbour);
                            float3 neighbourVelocity = Velocities[neighborIndex]; // Read current velocity
                            float3 velocity = Velocities[index]; // Read current velocity
                            viscosityForce += (neighbourVelocity - velocity) * SmoothingKernelPoly6(dst, ParticleSettings.SmoothingLength);

                        } while (SpatialGrid.TryGetNextValue(out neighborIndex, ref iterator));
                    }
                }

                // Write updated velocity to UpdatedVelocities
                UpdatedVelocities[index] = Velocities[index] + viscosityForce * ParticleSettings.ViscosityCoefficient * DeltaTime;
            }
        }

        [BurstCompile]
        public partial struct ResolveCollisionsApplyForcesJob : IJobEntity
        {
            [ReadOnly] public NativeArray<LocalToWorld> ParticlePositions;
            [ReadOnly] public NativeArray<Obstacle> Obstacles;

            [ReadOnly] public float DeltaTime;
            [ReadOnly] public Particle ParticleSettings;
            [ReadOnly] public float3 BoxSize;

            public NativeArray<float3> Velocities;

            public void Execute([EntityIndexInQuery] int index, ref LocalToWorld localToWorld)
            {
                float3 velocity = Velocities[index];

                float3 position = ParticlePositions[index].Position;
                position += Velocities[index] * DeltaTime;

                float3 halfSize = BoxSize / 2;
                float3 edgeDst = halfSize - math.abs(position);

                if (edgeDst.x <= 0)
                {
                    position.x = halfSize.x * math.sign(position.x);
                    velocity.x *= -ParticleSettings.CollisionDamping;
                }
                if (edgeDst.y <= 0)
                {
                    position.y = halfSize.y * math.sign(position.y);
                    velocity.y *= -ParticleSettings.CollisionDamping;
                }
                if (edgeDst.z <= 0)
                {
                    position.z = halfSize.z * math.sign(position.z);
                    velocity.z *= -ParticleSettings.CollisionDamping;
                }

                foreach (var obstacle in Obstacles)
                {
                    quaternion inverseRotation = math.inverse(obstacle.Rotation);
                    float3 localPosition = math.mul(inverseRotation, position - obstacle.Center);
                    float3 localVelocity = math.mul(inverseRotation, velocity);

                    float3 obstacleHalfSize = obstacle.Size * 0.5f + 0.5f;
                    float3 obstacleEdgeDst = obstacleHalfSize - math.abs(localPosition);

                    if (obstacleEdgeDst.x >= 0 && obstacleEdgeDst.y >= 0 && obstacleEdgeDst.z >= 0)
                    {
                        if (obstacleEdgeDst.x < obstacleEdgeDst.y && obstacleEdgeDst.x < obstacleEdgeDst.z)
                        {
                            localPosition.x = obstacleHalfSize.x * math.sign(localPosition.x);
                            localVelocity.x *= -ParticleSettings.CollisionDamping;
                        }
                        else if (obstacleEdgeDst.y < obstacleEdgeDst.z)
                        {
                            localPosition.y = obstacleHalfSize.y * math.sign(localPosition.y);
                            localVelocity.y *= -ParticleSettings.CollisionDamping;
                        }
                        else
                        {
                            localPosition.z = obstacleHalfSize.z * math.sign(localPosition.z);
                            localVelocity.z *= -ParticleSettings.CollisionDamping;
                        }

                        position = math.mul(obstacle.Rotation, localPosition) + obstacle.Center;
                        velocity = math.mul(obstacle.Rotation, localVelocity);
                    }
                }

                Velocities[index] = velocity;

                localToWorld = new LocalToWorld
                {
                    Value = float4x4.TRS(position, quaternion.identity, new float3(1, 1, 1))
                };
            }
        }

        private static float SmoothingKernelPoly6(float dst, float radius)
        {
            if (dst < radius)
            {
                float scale = 315 / (64 * math.PI * math.pow(math.abs(radius), 9));
                float v = radius * radius - dst * dst;
                return v * v * v * scale;
            }
            return 0;
        }

        private static float SpikyKernelPow2(float dst, float radius)
        {
            if (dst < radius)
            {
                float scale = 15 / (2 * math.PI * math.pow(radius, 5));
                float v = radius - dst;
                return v * v * scale;
            }
            return 0;
        }

        private static float SpikyKernelPow3(float dst, float radius)
        {
            if (dst < radius)
            {
                float scale = 15 / (math.PI * math.pow(radius, 6));
                float v = radius - dst;
                return v * v * v * scale;
            }
            return 0;
        }

        private static float DerivativeSpikyPow2(float dst, float radius)
        {
            if (dst <= radius)
            {
                float scale = 15 / (math.pow(radius, 5) * math.PI);
                float v = radius - dst;
                return -v * scale;
            }
            return 0;
        }

        private static float DerivativeSpikyPow3(float dst, float radius)
        {
            if (dst <= radius)
            {
                float scale = 45 / (math.pow(radius, 6) * math.PI);
                float v = radius - dst;
                return -v * v * scale;
            }
            return 0;
        }

        private static int3 SpatialHash(float3 position, float cellSize)
        {
            return new int3(
                (int)math.floor(position.x / cellSize),
                (int)math.floor(position.y / cellSize),
                (int)math.floor(position.z / cellSize)
            );
        }

    }
}
