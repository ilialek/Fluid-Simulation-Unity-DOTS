using Unity.Entities;
using Unity.Transforms;
using Unity.Mathematics;

namespace SPH
{
    [UpdateInGroup(typeof(SimulationSystemGroup))]
    public partial struct TransformAndCollectObstaclesSystem : ISystem
    {
        public void OnUpdate(ref SystemState state)
        {
            // Query entities with Obstacle component and dynamic buffer
            foreach (var (obstacle, localTransform) in SystemAPI
                     .Query<RefRW<Obstacle>, RefRW<LocalTransform>>()
                     )
            {
                // Rotate the obstacle using its RotateSpeed
                localTransform.ValueRW = localTransform.ValueRO.RotateY(obstacle.ValueRW.RotateSpeed * SystemAPI.Time.DeltaTime);

                // Update the obstacle's rotation in the Obstacle component
                obstacle.ValueRW.Rotation = localTransform.ValueRO.Rotation;

                // Update the dynamic buffer for this obstacle
                //if (SystemAPI.HasBuffer<ObstacleBufferElement>(entity))
                //{
                //    var buffer = SystemAPI.GetBuffer<ObstacleBufferElement>(entity);

                //    // Clear and update the buffer with new obstacle data
                //    buffer.Clear();
                //    buffer.Add(new ObstacleBufferElement
                //    {
                //        Size = obstacle.ValueRO.Size,
                //        Center = localTransform.ValueRO.Position,
                //        Rotation = localTransform.ValueRO.Rotation
                //    });
                //}
            }
        }
    }
}
