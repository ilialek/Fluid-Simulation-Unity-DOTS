using Unity.Entities;
using Unity.Mathematics;

namespace SPH
{
    [InternalBufferCapacity(4)]
    public struct ObstacleBufferElement : IBufferElementData
    {
        public float3 Size;
        public float3 Center;
        public quaternion Rotation;
    }
}
