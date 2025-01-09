using Unity.Burst;
using Unity.Collections;
using Unity.Collections.LowLevel.Unsafe;
using Unity.Entities;
using Unity.Jobs;
using Unity.Mathematics;
using Unity.Transforms;

namespace SPH
{
    [RequireMatchingQueriesForUpdate]
    [BurstCompile]
    public partial struct ParticleSpawnSystem : ISystem
    {
        private EntityQuery _particleQuery;
        public void OnCreate(ref SystemState state)
        {
            _particleQuery = new EntityQueryBuilder(Allocator.Temp).WithAll<Particle, LocalTransform>().Build(ref state);
        }

        public void OnUpdate(ref SystemState state)
        {
            ComponentLookup<LocalToWorld> localToWorldLookup = SystemAPI.GetComponentLookup<LocalToWorld>();
            EntityCommandBuffer ecb = new EntityCommandBuffer(Allocator.Temp);
            WorldUnmanaged world = state.World.Unmanaged;

            foreach (var (particleSpawn, particleSpawnLocalToWorld, entity) in
                     SystemAPI.Query<RefRO<ParticleSpawn>, RefRO<LocalToWorld>>()
                         .WithEntityAccess())
            {
                NativeArray<Entity> particleEntities = CollectionHelper.CreateNativeArray<Entity, RewindableAllocator>(particleSpawn.ValueRO.Count, ref world.UpdateAllocator);

                state.EntityManager.Instantiate(particleSpawn.ValueRO.ParticlePrefab, particleEntities);

                //float3 Origin = particleSpawnLocalToWorld.ValueRO.Position;
                float3 Origin = new float3(particleSpawnLocalToWorld.ValueRO.Position.x, particleSpawnLocalToWorld.ValueRO.Position.y + particleSpawn.ValueRO.BoxSize.y / 4, particleSpawnLocalToWorld.ValueRO.Position.z);
                float3 BoxSize = particleSpawn.ValueRO.BoxSize;

                SetParticleLocalToWorld setParticleLocalToWorldJob = new SetParticleLocalToWorld
                {
                    LocalToWorldFromEntity = localToWorldLookup,
                    Entities = particleEntities,
                    Origin = Origin,
                    BoxSize = BoxSize
                };
                state.Dependency = setParticleLocalToWorldJob.Schedule(particleSpawn.ValueRO.Count, 64, state.Dependency);
                state.Dependency.Complete();

                ecb.DestroyEntity(entity);
            }

            ecb.Playback(state.EntityManager);
            // All Prefabs are currently forced to TransformUsageFlags.Dynamic by default, which means boids get a LocalTransform
            // they don't need. As a workaround, remove the component at spawn-time.
            state.EntityManager.RemoveComponent<LocalTransform>(_particleQuery);
        }
    }

    [BurstCompile]
    struct SetParticleLocalToWorld : IJobParallelFor
    {
        [NativeDisableContainerSafetyRestriction]
        [NativeDisableParallelForRestriction]
        public ComponentLookup<LocalToWorld> LocalToWorldFromEntity;

        public NativeArray<Entity> Entities;
        public float3 Origin;

        public float3 BoxSize;

        public void Execute(int i)
        {
            Entity entity = Entities[i];

            Random random = new Random(((uint)(entity.Index + i + 1) * 0x9F6ABC1));

            float3 pos = Origin + random.NextFloat3(-BoxSize / 4, BoxSize / 4);

            //float3 pos = new float3(Origin.x + random.NextFloat(-BoxSize.x / 3, BoxSize.x / 3), Origin.y + random.NextFloat(-BoxSize.y / 3, BoxSize.y / 3), 0);

            LocalToWorld localToWorld = new LocalToWorld
            {
                Value = float4x4.TRS(pos, quaternion.identity, new float3(1, 1, 1))
            };

            // Assign the calculated LocalToWorld to the entity
            LocalToWorldFromEntity[entity] = localToWorld;
        }
    }

}
