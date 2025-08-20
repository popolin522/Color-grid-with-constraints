import numpy as np
import matplotlib.pyplot as plt

def create_facet_grid_pattern(vertex_gd=0.5, edge_gd=0.5, face_gd=0.5, blue_ratio=0.3):
    """
    Create a 7×7 grid pattern for a specific facet with uniform parameters across all facets.
    
    Parameters:
    - vertex_gd: Grafting density for vertices/corners (0.0-1.0)
    - edge_gd: Grafting density for edges (0.0-1.0)
    - face_gd: Grafting density for face/interior cells (0.0-1.0)
    - blue_ratio: Proportion of molecules that should be hydrophobic, or blue in the coloring problem (0.0-1.0)
    
    Returns:
    - grid: 7×7 numpy array with values 0 (white/removed), 1 (green/non-hydrophobic), 2 (blue/hydrophobic)
    """
    # Make molecules and assign them sequentially
    facet_molecules = list(range(49))
    grid_molecules = np.zeros((7, 7), dtype=int)
    
    # Arrange molecules in the grid (row by row, with IDs increasing left to right)
    for i in range(7):
        for j in range(7):
            index = i * 7 + j
            grid_molecules[i, j] = facet_molecules[index]
    
    # Define geometric regions based on the create_grid_pattern function
    corners = [(0,0), (0,1), (1,0), (0,6), (0,5), (1,6), (6,0), (5,0), (6,1), (6,6), (6,5), (5,6)]
    
    edges = []
    # Outer perimeter (excluding corners)
    for i in range(7):
        for j in range(7):
            if i == 0 or i == 6 or j == 0 or j == 6:
                if (i,j) not in corners:
                    edges.append((i,j))
    
    # Inner perimeter (excluding corners)
    for i in range(1,6):
        for j in range(1,6):
            if i == 1 or i == 5 or j == 1 or j == 5:
                if (i,j) not in corners:
                    edges.append((i,j))
    
    # Face cells are the remaining interior cells
    faces = [(i,j) for i in range(2,5) for j in range(2,5)]
    
    # Calculate the total number of chains and color ratios based on input parameters
    vertex_chains = int(vertex_gd * len(corners))
    edge_chains = int(edge_gd * len(edges))
    face_chains = int(face_gd * len(faces))
    
    # Calculate total chains and hydrophobic counts
    total_chains = vertex_chains + edge_chains + face_chains
    

    hydrophobic_count = int(total_chains * blue_ratio)

    # Calculate non-hydrophobic count
    non_hydrophobic_count = total_chains - hydrophobic_count
    
    # Create the grid pattern
    grid = np.zeros((7, 7), dtype=int)
    
    # Initialize sets to track molecules
    
    # Apply blue (hydrophobic) allocation with priority: corners, edges, faces
    blue_count = 0
    
    # Corner groups for symmetry
    corner_groups = [
        [(0,0), (6,6)],  # Diagonal corners group 1
        [(0,6), (6,0)],  # Diagonal corners group 2
        [(0,1), (6,5)],  # Adjacent to corners with symmetry
        [(1,0), (5,6)],
        [(0,5), (6,1)],
        [(1,6), (5,0)]
    ]
    
    # Keep track of number of colored sites
    current_vertex_count = 0
    current_edge_count = 0
    current_face_count = 0
    
    # 1. First priority: Fill corners with blue
    for group in corner_groups:
        for pos in group:
            if current_vertex_count >= vertex_chains or blue_count >= hydrophobic_count:
                break
            
            grid[pos] = 2  # Blue
            blue_count += 1
            current_vertex_count += 1
        if current_vertex_count >= vertex_chains:
            break
    
    # 2. Second priority: Fill edges with blue
    edge_groups = [
        [(6,3), (3,6), (0,3), (3,0)],  # Middle of outer edges
        [(1,1), (5,5), (1,5), (5,1)],
        [(5,2), (1,4)],
        [(0,2), (0,4), (2,0), (4,0)], [(6,2), (6,4), (2,6), (4,6)],  # Other outer edges
        [(1,2), (5,4), (1,3), (5,3)],  # Inner horizontal edges
        [(2,1), (4,5), (4,1)], [(2,5), (3,1), (3,5)]   # Inner vertical edges
    ]
    
    sites_at_high_curvature_region = vertex_chains + edge_chains
    
    for group in edge_groups:
        available_positions = [pos for pos in group if pos in edges and grid[pos] == 0]
        positions_to_fill = min(len(available_positions), hydrophobic_count - blue_count)

        if positions_to_fill <= 0:
            break
            
        # Prioritize maintaining symmetry
        if positions_to_fill >= 2 and positions_to_fill % 2 == 0 and len(available_positions) % 2 == 0:
            # Fill in pairs to maintain symmetry
            for i in range(0, positions_to_fill, 2):
                if i+1 < len(available_positions):
                    pos1, pos2 = available_positions[i], available_positions[i+1]
                    grid[pos1] = 2
                    grid[pos2] = 2
                    blue_count += 2
                    current_edge_count += 2
        else:
            # Fill as many as possible
            for i in range(positions_to_fill):
                if i < len(available_positions):
                    pos = available_positions[i]
                    grid[pos] = 2
                    blue_count += 1
                    current_edge_count += 1
        if current_edge_count >= edge_chains or blue_count >= hydrophobic_count:
            break
    
    # 3. Third priority: Fill faces with blue in a symmetrical pattern
    face_groups = [
        [(3,3)],  # Center
        [(2,2), (4,4), (4,2), (2,4)],  # Diagonal from center
        [(2,3), (4,3), (3,2), (3,4)]   # Adjacent to center
    ]
    
    total_chains_target = sites_at_high_curvature_region + face_chains
    remaining_blue = min(total_chains_target - blue_count, hydrophobic_count - blue_count)
    
    for group in face_groups:
        available_positions = [pos for pos in group if pos in faces and grid[pos] == 0]
        
        if len(available_positions) <= remaining_blue:
            # Fill the whole group
            for pos in available_positions:
                grid[pos] = 2  # Blue
                blue_count += 1
                current_face_count += 1
                remaining_blue -= 1
        else:
            # Not enough remaining spots, but aim for symmetry
            if len(group) == 1:  # Single center cell
                if remaining_blue >= 1:
                    grid[group[0]] = 2
                    blue_count += 1
                    current_face_count += 1
                    remaining_blue -= 1
            elif remaining_blue >= 2:  # Fill pairs when possible
                pairs_to_fill = min(remaining_blue // 2, len(available_positions) // 2)
                for i in range(pairs_to_fill):
                    if 2*i+1 < len(available_positions):
                        pos1, pos2 = available_positions[2*i], available_positions[2*i+1]
                        grid[pos1] = 2
                        grid[pos2] = 2
                        blue_count += 2
                        current_face_count += 2
                        remaining_blue -= 2
            
        if remaining_blue <= 0:
            break
        if current_face_count >= face_chains:
            break
    
    # Now distribute green (non-hydrophobic) with similar symmetry patterns
    green_count = 0
    target_green_remaining = non_hydrophobic_count
    
    # First try to place green in face positions that aren't already blue
    for group in face_groups:
        available_pos = [pos for pos in group if grid[pos] == 0]
        positions_to_fill = min(len(available_pos), target_green_remaining)
        
        if positions_to_fill > 0:
            if positions_to_fill == len(available_pos):
                # Fill all available positions
                for pos in available_pos:
                    grid[pos] = 1  # Green
                    green_count += 1
                    current_face_count += 1
                    target_green_remaining -= 1
            elif positions_to_fill % 2 == 0 and len(available_pos) % 2 == 0:
                # Fill in pairs to maintain symmetry
                for i in range(0, positions_to_fill, 2):
                    if i+1 < len(available_pos):
                        pos1, pos2 = available_pos[i], available_pos[i+1]
                        grid[pos1] = 1
                        grid[pos2] = 1
                        green_count += 2
                        current_face_count += 2
                        target_green_remaining -= 2
            else:
                # Fill as many as possible
                for i in range(positions_to_fill):
                    pos = available_pos[i]
                    grid[pos] = 1
                    green_count += 1
                    current_face_count += 1
                    target_green_remaining -= 1
        if current_face_count >= face_chains:
            break
    
    # Then corner positions
    for group in corner_groups:
        available_pos = [pos for pos in group if grid[pos] == 0]
        positions_to_fill = min(len(available_pos), target_green_remaining)
        
        if positions_to_fill > 0:
            if positions_to_fill == len(available_pos):
                # Fill all available positions
                for pos in available_pos:
                    grid[pos] = 1  # Green
                    green_count += 1
                    current_vertex_count += 1
                    target_green_remaining -= 1
            elif positions_to_fill % 2 == 0 and len(available_pos) % 2 == 0:
                # Fill in pairs to maintain symmetry
                for i in range(0, positions_to_fill, 2):
                    if i+1 < len(available_pos):
                        pos1, pos2 = available_pos[i], available_pos[i+1]
                        grid[pos1] = 1
                        grid[pos2] = 1
                        green_count += 2
                        current_vertex_count += 2
                        target_green_remaining -= 2
            else:
                # Fill as many as possible
                for i in range(positions_to_fill):
                    pos = available_pos[i]
                    grid[pos] = 1
                    green_count += 1
                    current_vertex_count += 1
                    target_green_remaining -= 1
        
        if target_green_remaining <= 0:
            break
        if current_vertex_count >= vertex_chains:
            break

    # Finally edge positions
    for group in reversed(edge_groups):
        available_pos = [pos for pos in group if grid[pos] == 0]
        positions_to_fill = min(len(available_pos), target_green_remaining)
        
        if positions_to_fill > 0:
            if positions_to_fill == len(available_pos):
                # Fill all available positions
                for pos in available_pos:
                    grid[pos] = 1  # Green
                    green_count += 1
                    current_edge_count += 1
                    target_green_remaining -= 1
            elif positions_to_fill % 2 == 0 and len(available_pos) % 2 == 0:
                # Fill in pairs to maintain symmetry
                for i in range(0, positions_to_fill, 2):
                    if i+1 < len(available_pos):
                        pos1, pos2 = available_pos[i], available_pos[i+1]
                        grid[pos1] = 1
                        grid[pos2] = 1
                        green_count += 2
                        current_edge_count += 2
                        target_green_remaining -= 2
            else:
                # Fill as many as possible
                for i in range(positions_to_fill):
                    pos = available_pos[i]
                    grid[pos] = 1
                    green_count += 1
                    current_edge_count += 1
                    target_green_remaining -= 1
        
        if target_green_remaining <= 0:
            break
        if current_edge_count >= edge_chains:
            break
    
    return grid

def create_sphere_grid_visualization(grid, title="", show_stats=True):
    """
    Create a sphere-based visualization of the 7x7 grid pattern.
    
    Parameters:
    - grid: 7x7 numpy array with values 0 (white), 1 (green), 2 (blue)
    - ax: matplotlib axis to plot on
    - title: title for the subplot
    - show_stats: whether to show statistics text
    """
    # Define colors for spheres
    colors = ['grey', 'green', 'blue']
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(6, 6))
    # Create sphere positions in 7x7 grid
    sphere_radius = 0.45
    spacing = 1.0
    
    # Plot spheres for each position
    for i in range(7):
        for j in range(7):
            x = j * spacing
            y = (6-i) * spacing   #flip y-axis to match grid orientation
            
            # Get color based on grid value
            color = colors[grid[i, j]]
            
            # Create circle (sphere in 2D)
            circle = plt.Circle((x, y), sphere_radius, color=color, edgecolor='black', linewidth=1)
            ax.add_patch(circle)
    
    # Set axis limits and properties
    ax.set_xlim(-0.5, 6.5)
    ax.set_ylim(-0.5, 6.5)
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Add title
    if title:
        ax.set_title(title, fontsize=10)#, pad=10)
    
    # Add statistics if requested
    if show_stats:
        # Calculate statistics
        white_count = np.sum(grid == 0)
        green_count = np.sum(grid == 1)
        blue_count = np.sum(grid == 2)
        total_grafted = green_count + blue_count
        
        # Calculate actual hydrophobic ratio
        actual_h_ratio = blue_count / total_grafted if total_grafted > 0 else 0
        
        # Calculate total grafting density
        total_gd = total_grafted / 49.0
        
        # Add text annotations
        # ax.text(6.2, 6.2, f'GD: {total_gd:.2f}', fontsize=9, ha='right', va='top', 
        #         bbox=dict(boxstyle="round,pad=0.9", facecolor='white', alpha=0.8))
        x_center = 0.5 * sum(ax.get_xlim())

        # ax.text(x_center, -0.5, f'$\Gamma$,H = ({total_gd:.2f},{actual_h_ratio:.2f})', fontsize=7, ha='center', va='bottom',
        #         bbox=dict(boxstyle="square,pad=0.1", facecolor='white', alpha=0.8))
    plt.tight_layout()
    return fig


def create_sample_visualization(vertex_gd, edge_gd, face_gd, blue_ratio):
    """
    Create a sample visualization to demonstrate the function.
    """
    
    # Create grid pattern
    grid = create_facet_grid_pattern(
        vertex_gd, edge_gd, face_gd, blue_ratio
    )
    
    # Create visualization
    title = f"Design constraint: v={vertex_gd}, e={edge_gd}, f={face_gd}, H={blue_ratio}"
    fig = create_sphere_grid_visualization(grid, title=title, show_stats=True)
    
    # Save the figure
    fig.savefig('sample_color_grid.png', dpi=300, bbox_inches='tight')
    
    plt.show()
    return fig


if __name__ == "__main__":
    #design parameters. change here
    vertex_gd = 1.0
    edge_gd = 0.5
    face_gd = 0.5
    blue_ratio = 0.5
    create_sample_visualization(vertex_gd, edge_gd, face_gd, blue_ratio)