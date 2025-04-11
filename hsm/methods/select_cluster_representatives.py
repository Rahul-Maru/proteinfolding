import py4cytoscape as p4c
import pandas as pd
import os
import re
import sys
from collections import defaultdict, deque

# --- Configuration ---
# IMPORTANT: Adjust these values to match your setup

# !! STEP 1: Define the attribute name we WANT to create for components !!
# We will explicitly tell the command to use this name.
CLUSTER_ID_ATTRIBUTE = 'MyComponentID' # <-- Let's try creating our own named column

# Attribute to store the result (True for representatives, False otherwise)
REPRESENTATIVE_ATTRIBUTE = 'is_cluster_representative'
# Set to True to also select the non-representative nodes in Cytoscape UI
SELECT_REPRESENTATIVES_IN_UI = True

# --- Paths ---
# ... (existing path setup code) ...
try:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    hsm_dir = os.path.dirname(script_dir) # Assumes script is in hsm/methods
    PDB_DIR = os.path.join(hsm_dir, 'pdbs')
    BCHAIN_DIR = os.path.join(hsm_dir, 'bchains')
except NameError:
    print("Warning: Could not determine script location automatically...")
    PDB_DIR = 'hsm/pdbs'
    BCHAIN_DIR = 'hsm/bchains'

# Check if data directories exist
if not os.path.isdir(PDB_DIR):
    print(f"Error: PDB directory not found: {PDB_DIR}")
    sys.exit(1)
if not os.path.isdir(BCHAIN_DIR):
    print(f"Error: BChain directory not found: {BCHAIN_DIR}")
    sys.exit(1)

# --- Helper Functions ---
# ... (get_pdb_resolution and get_bchain_filesize functions remain the same) ...
def get_pdb_resolution(pdb_id):
    """
    Parses a PDB file in PDB_DIR to extract resolution from REMARK 2.
    Returns resolution as float (lower is better), or infinity if not found/error.
    """
    pdb_file_path = os.path.join(PDB_DIR, f"{pdb_id}.pdb")
    try:
        with open(pdb_file_path, 'r') as f:
            for line in f:
                if line.startswith("REMARK   2 RESOLUTION."):
                    # Extract floating point number
                    match = re.search(r'(\d+\.\d+)\s+ANGSTROMS', line)
                    if match:
                        resolution = float(match.group(1))
                        # print(f" Found resolution {resolution} for {pdb_id}")
                        return resolution
                    # Handle cases like "NOT APPLICABLE" or "NULL" explicitly if needed
                    if "NOT APPLICABLE" in line or "NULL" in line:
                         # print(f" Resolution not applicable for {pdb_id}")
                         return float('inf') # Assign infinite resolution (worst)
            # print(f" Resolution REMARK not found for {pdb_id}")
            return float('inf') # Return worst resolution if REMARK 2 not found or no value
    except FileNotFoundError:
        # print(f" Warning: PDB file not found: {pdb_file_path}")
        return float('inf') # Worst resolution if file not found
    except Exception as e:
        print(f" Warning: Error reading/parsing PDB {pdb_file_path}: {e}")
        return float('inf') # Worst resolution on error

def get_bchain_filesize(node_name):
    """
    Gets the file size for a bchain file in BCHAIN_DIR.
    Assumes filename is {node_name}.pdb. Returns size in bytes or -1 if not found/error.
    """
    # Adjust extension if necessary (e.g., .ent, .chain)
    bchain_file_path = os.path.join(BCHAIN_DIR, f"{node_name}.pdb")
    try:
        if os.path.exists(bchain_file_path):
            size = os.path.getsize(bchain_file_path)
            # print(f" Found bchain size {size} for {node_name}")
            return size
        else:
            # print(f" Warning: BChain file not found: {bchain_file_path}")
            return -1 # Return -1 (worst size) if file not found
    except Exception as e:
        print(f" Warning: Error getting size for {bchain_file_path}: {e}")
        return -1 # Return -1 on error

def find_connected_components(nodes, edges):
    """
    Identifies connected components in a network using BFS.

    Args:
        nodes: Dictionary mapping suid to node info (name, degree)
        edges: List of tuples (source, target)

    Returns:
        Dictionary mapping component_id to list of node dictionaries in that component
    """
    # Create adjacency list
    adjacency = defaultdict(list)
    for source, target in edges:
        adjacency[source].append(target)
        adjacency[target].append(source) # Bidirectional

    # Track visited nodes
    visited = set()
    components = defaultdict(list)
    component_id = 0

    # Run BFS from each unvisited node
    for suid in nodes.keys():
        if suid in visited:
            continue

        # Start a new component with BFS
        component = []
        queue = deque([suid])
        visited.add(suid)

        # BFS to find all connected nodes
        while queue:
            current = queue.popleft()
            component.append(current)

            # Visit all unvisited neighbors
            for neighbor in adjacency[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)

        # Store this component's nodes (convert SUIDs to node objects)
        components[component_id] = [
            {'suid': suid, 'name': nodes[suid]['name'], 'degree': nodes[suid]['degree']}
            for suid in component if suid in nodes
        ]
        component_id += 1

    # Check for isolated nodes (nodes with no edges)
    for suid, node_info in nodes.items():
        if suid not in visited:  # This node had no edges
            components[component_id] = [{'suid': suid, 'name': node_info['name'], 'degree': node_info['degree']}]
            component_id += 1
            visited.add(suid)

    return components

def select_representatives(network=None):
    """
    Identifies connected components manually, selects representative nodes,
    and selects non-representative nodes in the Cytoscape UI.

    Args:
        network (str, int, None): Network name, SUID, or None for current network.
    """
    try:
        # --- Connection and Network Setup ---
        print("Connecting to Cytoscape...")
        p4c.cytoscape_ping()
        print("Connection successful.")

        network_suid = p4c.get_network_suid(network)
        network_name = p4c.get_network_name(network_suid)
        print(f"Processing network: {network_name} (SUID: {network_suid})")

        # --- Fetch Nodes Using Alternative Methods ---
        print("Fetching node data...")

        # Get all nodes
        all_nodes = p4c.get_all_nodes(network=network_suid)
        print(f"Found {len(all_nodes)} nodes.")

        # Get node names (could also use get_node_attribute)
        try:
            node_names = p4c.get_node_names(network=network_suid)
        except:
            print("Warning: Could not get node names using get_node_names, trying alternative...")
            # Alternative: get the shared name or name attribute directly
            node_names = p4c.get_node_attribute(names=all_nodes, column='name', network=network_suid)

        # Get node degree (adjacency count)
        try:
            # Try a direct approach to get degree
            node_degrees = p4c.get_node_attribute(names=all_nodes, column='degree', network=network_suid)
        except:
            print("Warning: Could not get node degrees, calculating from edges...")
            # If no degree column exists, calculate from edges
            # Initialize with zero degrees
            node_degrees = {suid: 0 for suid in all_nodes}
            # We'll add to this when we process edges

        # Build node info dictionary
        nodes = {}
        for i, suid in enumerate(all_nodes):
            nodes[suid] = {
                'name': node_names.get(suid, f"Node_{suid}"),
                'degree': node_degrees.get(suid, 0)
            }

        # --- Fetch Edges ---
        print("Fetching edge data...")
        # Get all edges
        all_edges = p4c.get_all_edges(network=network_suid)
        print(f"Found {len(all_edges)} edges.")

        # Get edge source and target
        edge_source = p4c.get_edge_property(edge_names=all_edges, property_name="source", network=network_suid)
        edge_target = p4c.get_edge_property(edge_names=all_edges, property_name="target", network=network_suid)

        # Create list of edges as (source, target) tuples
        edges = []
        for edge in all_edges:
            if edge in edge_source and edge in edge_target:
                source = edge_source[edge]
                target = edge_target[edge]
                edges.append((source, target))

                # If we couldn't get degree attribute earlier, calculate it now
                if 'degree' not in node_degrees:
                    node_degrees[source] = node_degrees.get(source, 0) + 1
                    node_degrees[target] = node_degrees.get(target, 0) + 1
                    # Also update our nodes dictionary
                    if source in nodes:
                        nodes[source]['degree'] = node_degrees[source]
                    if target in nodes:
                        nodes[target]['degree'] = node_degrees[target]

        # --- Find Connected Components Manually ---
        print("Finding connected components based on edge data...")
        components = find_connected_components(nodes, edges)
        print(f"Found {len(components)} connected components.")

        # Store all node SUIDs for selection logic later
        all_node_suids = set(all_nodes)

        # --- Select Representatives ---
        print("Selecting representative for each component...")
        representatives = {}  # {component_id: representative_suid}

        for component_id, component_nodes in components.items():
            if not component_nodes:
                continue

            # Gather data for each node in this component
            cluster_candidates = []
            for node in component_nodes:
                node_name = node['name']
                pdb_id = node_name.split('_')[0] if '_' in node_name else node_name
                resolution = get_pdb_resolution(pdb_id)
                filesize = get_bchain_filesize(node_name)

                cluster_candidates.append({
                   'suid': node['suid'],
                   'name': node_name,
                   'degree': node['degree'],
                   'resolution': resolution,
                   'filesize': filesize
                })

            # Sort candidates based on criteria
            sorted_candidates = sorted(
                cluster_candidates,
                key=lambda x: (
                    x['degree'],            # Maximize degree
                    -x['resolution'],       # Minimize resolution (so negate)
                    x['filesize'],          # Maximize filesize
                    x['name']               # Alphabetical tie-breaker
                ),
                reverse=True
            )

            # Select best candidate as representative
            if sorted_candidates:
                representative_node = sorted_candidates[0]
                representatives[component_id] = representative_node['suid']
                # Debug print for first few components
                if component_id < 5:  # Show details for the first 5 components
                    print(f"  Component {component_id} ({len(component_nodes)} nodes)")
                    print(f"    Representative: {representative_node['name']} (Degree: {representative_node['degree']},",
                          f"Resolution: {representative_node['resolution']},",
                          f"File Size: {representative_node['filesize']})")

            # Print progress for large networks
            if component_id > 0 and component_id % 100 == 0:
                print(f"  Processed {component_id}/{len(components)} components...")

        # Gather all representative SUIDs
        representative_suids = set(representatives.values())
        print(f"Total representative nodes identified: {len(representative_suids)}")

        # --- Update Node Table ---
        print(f"Updating node table with '{REPRESENTATIVE_ATTRIBUTE}' attribute...")
        # Create a dictionary mapping node SUID to whether it's a representative
        is_representative = {suid: (suid in representative_suids) for suid in all_nodes}

        # Use set_node_attribute to update the Cytoscape table
        p4c.set_node_attribute(node_names=all_nodes, column=REPRESENTATIVE_ATTRIBUTE,
                              values=is_representative, network=network_suid)
        print("Node table updated.")

        # --- Select Non-Representatives in UI ---
        if SELECT_REPRESENTATIVES_IN_UI:
            non_representative_suids = all_node_suids - representative_suids
            if non_representative_suids:
                print(f"Selecting {len(non_representative_suids)} non-representative nodes in Cytoscape UI...")
                p4c.clear_selection(type='nodes', network=network_suid)
                p4c.select_nodes(nodes=list(non_representative_suids), network=network_suid)
                print("Non-representative nodes selected.")
            else:
                print("No non-representative nodes found to select.")

        print("\nScript finished successfully.")
        return True

    except ImportError as e:
         print(f"Error: Missing required library. {e}")
         print("Please install py4cytoscape: pip install py4cytoscape")
         return False
    except Exception as e:
        print(f"\nAn error occurred: {e}")
        print("Please ensure Cytoscape is running, the network is loaded,")
        print("and the specified paths and attributes are correct.")
        import traceback
        traceback.print_exc() # Print detailed traceback for debugging
        return False

def main():
    """Main function to run the cluster representative selection process."""
    # Check if data directories exist
    if not os.path.isdir(PDB_DIR):
        print(f"Error: PDB directory not found: {PDB_DIR}")
        return 1
    if not os.path.isdir(BCHAIN_DIR):
        print(f"Error: BChain directory not found: {BCHAIN_DIR}")
        return 1

    # Get target network from command line if specified
    target_network = sys.argv[1] if len(sys.argv) > 1 else None

    # Run the representative selection process
    success = select_representatives(network=target_network)

    # Return appropriate exit code
    return 0 if success else 1

# --- Run the script ---
if __name__ == '__main__':
    sys.exit(main())