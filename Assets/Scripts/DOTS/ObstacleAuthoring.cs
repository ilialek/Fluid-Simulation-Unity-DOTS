using Unity.Entities;
using UnityEngine;
using Unity.Mathematics;

namespace SPH
{
    public class ObstacleAuthoring : MonoBehaviour
    {
        public float RotateSpeed;

        class Baker : Baker<ObstacleAuthoring>
        {
            public override void Bake(ObstacleAuthoring authoring)
            {
                var entity = GetEntity(TransformUsageFlags.Dynamic);
                AddComponent(entity, new Obstacle
                {
                    RotateSpeed = authoring.RotateSpeed,
                    Size = authoring.transform.localScale,
                    Center = authoring.transform.position,
                    Rotation = authoring.transform.rotation
                });
            }
        }
    }

    public struct Obstacle : IComponentData
    {
        public float RotateSpeed;
        public float3 Size;
        public float3 Center;
        public quaternion Rotation;
    }
}
