using UnityEngine;

public class WireframeBoxGizmos : MonoBehaviour
{
    public Vector3 boxSize = new Vector3(1, 1, 1); // Dimensions of the box
    public Color boxColor = Color.white;          // Color of the edges

    void OnDrawGizmos()
    {
        Gizmos.color = boxColor;

        // Calculate box corners (relative to the object's position)
        Vector3[] corners = {
            transform.position + new Vector3(-boxSize.x, -boxSize.y, -boxSize.z) / 2,
            transform.position + new Vector3(boxSize.x, -boxSize.y, -boxSize.z) / 2,
            transform.position + new Vector3(boxSize.x, boxSize.y, -boxSize.z) / 2,
            transform.position + new Vector3(-boxSize.x, boxSize.y, -boxSize.z) / 2,
            transform.position + new Vector3(-boxSize.x, -boxSize.y, boxSize.z) / 2,
            transform.position + new Vector3(boxSize.x, -boxSize.y, boxSize.z) / 2,
            transform.position + new Vector3(boxSize.x, boxSize.y, boxSize.z) / 2,
            transform.position + new Vector3(-boxSize.x, boxSize.y, boxSize.z) / 2
        };

        // Draw edges
        int[,] edges = {
            {0, 1}, {1, 2}, {2, 3}, {3, 0}, // Bottom face
            {4, 5}, {5, 6}, {6, 7}, {7, 4}, // Top face
            {0, 4}, {1, 5}, {2, 6}, {3, 7}  // Sides
        };

        for (int i = 0; i < edges.GetLength(0); i++)
        {
            Gizmos.DrawLine(corners[edges[i, 0]], corners[edges[i, 1]]);
        }
    }
}
