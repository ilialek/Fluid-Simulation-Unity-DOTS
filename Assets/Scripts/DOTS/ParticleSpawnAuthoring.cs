using Unity.Entities;
using UnityEngine;
using Unity.Mathematics;

namespace SPH
{
    public class ParticleSpawnAuthoring : MonoBehaviour
    {
        public GameObject ParticlePrefab;

        public int Count;

        public Vector3 BoxSize = new Vector3(1, 1, 1); // Dimensions of the box

        class Baker : Baker<ParticleSpawnAuthoring>
        {
            public override void Bake(ParticleSpawnAuthoring authoring)
            {
                var entity = GetEntity(TransformUsageFlags.Renderable);
                AddComponent(entity, new ParticleSpawn
                {
                    ParticlePrefab = GetEntity(authoring.ParticlePrefab, TransformUsageFlags.Renderable | TransformUsageFlags.WorldSpace),
                    Count = authoring.Count,
                    BoxSize = authoring.BoxSize
                });
            }
        }
    }

    public struct ParticleSpawn : IComponentData
    {
        public Entity ParticlePrefab;
        public int Count;
        public float3 BoxSize;
    }
}
