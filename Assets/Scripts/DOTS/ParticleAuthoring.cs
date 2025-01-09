using Unity.Entities;
using UnityEngine;
using Unity.Mathematics;

namespace SPH
{
    public class ParticleAuthoring : MonoBehaviour
    {
        [HideInInspector] public float Density;
        [HideInInspector] public float Pressure;
        [HideInInspector] public float3 PressureForce;
        [HideInInspector] public Vector3 Velocity = Vector3.zero;
        public float Mass;
        public float ReferenceDensity;
        public float Stiffness;
        public float SmoothingLength;
        public float CellRadius;
        public float CollisionDamping;
        public float Gravity;
        public float TimeScale;
        public float ViscosityCoefficient;
        public float NearPressureMultiplier;
        public int MaximumParticlesPerCell;

        class Baker : Baker<ParticleAuthoring>
        {
            public override void Bake(ParticleAuthoring authoring)
            {
                var entity = GetEntity(TransformUsageFlags.Renderable);
                AddComponent(entity, new Particle
                {
                    Density = authoring.Density,
                    Pressure = authoring.Pressure,
                    PressureForce = authoring.PressureForce,
                    Velocity = authoring.Velocity,
                    Mass = authoring.Mass,
                    ReferenceDensity = authoring.ReferenceDensity,
                    Stiffness = authoring.Stiffness,
                    SmoothingLength = authoring.SmoothingLength,
                    CellRadius = authoring.CellRadius,
                    CollisionDamping = authoring.CollisionDamping,
                    Gravity = authoring.Gravity,
                    TimeScale = authoring.TimeScale,
                    ViscosityCoefficient = authoring.ViscosityCoefficient,
                    NearPressureMultiplier = authoring.NearPressureMultiplier,
                    MaximumParticlesPerCell = authoring.MaximumParticlesPerCell
                });
            }
        }
    }

    public struct Particle : IComponentData
    {
        public float Density;
        public float Pressure;
        public float3 PressureForce;
        public float3 Velocity;
        public float Mass;
        public float ReferenceDensity;
        public float Stiffness;
        public float SmoothingLength;
        public float CellRadius;
        public float CollisionDamping;
        public float Gravity;
        public float TimeScale;
        public float ViscosityCoefficient;
        public float NearPressureMultiplier;
        public int MaximumParticlesPerCell;
    }
}
