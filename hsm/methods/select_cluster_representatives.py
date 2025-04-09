import py4cytoscape as p4c
import pandas as pd
import os
import re
import sys
from collections import defaultdict, deque
import time # Import time for potential delay

# --- Configuration ---
# IMPORTANT: Adjust these values to match your setup

# !! STEP 1: Define the attribute name we WANT to create for components !!
# We will explicitly tell the command to use this name.
CLUSTER_ID_ATTRIBUTE = 'MyComponentID' # <-- Let's try creating our own named column

# Attribute to store the result (True for representatives, False otherwise)
REPRESENTATIVE_ATTRIBUTE = 'is_cluster_representative'
# Set to True to also select the non-representative nodes in Cytoscape UI
SELECT_REPRESENTATIVES_IN_UI = True

# Try to determine paths relative to the script's location
try:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    hsm_dir = os.path.dirname(script_dir) # Assumes script is in hsm/methods
    PDB_DIR = os.path.join(hsm_dir, 'pdbs')
    BCHAIN_DIR = os.path.join(hsm_dir, 'bchains')
    # Add hsm_dir to sys.path if you need to import modules from there later
    # sys.path.insert(0, hsm_dir)
except NameError:
    # Fallback if __file__ is not defined (e.g., running interactively)
    print("Warning: Could not determine script location automatically. Using default relative paths.")
    print("Please ensure you run this script from the directory containing the 'hsm' folder,")
    print("or adjust PDB_DIR and BCHAIN_DIR manually.")
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
        nodes: List of node dictionaries with 'suid', 'name', 'degree'
        edges: List of edge dictionaries with 'source', 'target'
        
    Returns:
        Dictionary mapping component_id to list of nodes in that component
    """
    # Create adjacency list
    adjacency = defaultdict(list)
    for edge in edges:
        source = edge['source']
        target = edge['target']
        adjacency[source].append(target)
        adjacency[target].append(source) # Bidirectional
    
    # Track visited nodes
    visited = set()
    components = defaultdict(list)
    component_id = 0
    
    # Map SUID to node object for later
    suid_to_node = {node['suid']: node for node in nodes}
    
    # Run BFS from each unvisited node
    for node in nodes:
        suid = node['suid']
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
        
        # Store this component's nodes (convert SUIDs back to node objects)
        components[component_id] = [suid_to_node[suid] for suid in component if suid in suid_to_node]
        component_id += 1
    
    # Check for isolated nodes (nodes with no edges)
    for node in nodes:
        suid = node['suid']
        if suid not in visited:  # This node had no edges
            components[component_id] = [node]
            component_id += 1
            visited.add(suid)
    
    return components

# --- Main Logic ---

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

        # --- Fetch Nodes and Edges ---
        print("Fetching node data...")
        node_data_df = p4c.get_node_table(
            columns=['SUID', 'name', 'degree'],
            network=network_suid
        )
        # Convert to int and handle missing values
        node_data_df['SUID'] = node_data_df['SUID'].astype(int)
        node_data_df['degree'] = node_data_df['degree'].fillna(0)
        print(f"Found {len(node_data_df)} nodes.")
        
        print("Fetching edge data...")
        edge_data_df = p4c.get_edge_table(
            columns=['SUID', 'source', 'target'],
            network=network_suid
        )
        print(f"Found {len(edge_data_df)} edges.")
        
        # Convert DataFrames to lists of dictionaries for easier processing
        nodes = node_data_df.to_dict('records')
        edges = edge_data_df.to_dict('records')
        
        # --- Find Connected Components Manually ---
        print("Finding connected components based on edge data...")
        components = find_connected_components(nodes, edges)
        print(f"Found {len(components)} connected components.")
        
        # Store all node SUIDs for selection logic later
        all_node_suids = set(node_data_df['SUID'])
        
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
        update_df = pd.DataFrame({
            'SUID': list(node_data_df['SUID']),
            REPRESENTATIVE_ATTRIBUTE: [suid in representative_suids for suid in node_data_df['SUID']]
        })
        p4c.load_table_data(update_df, data_key_column='SUID', table='node', network=network_suid)
        print("Node table updated.")
        
        # --- Select Non-Representatives in UI ---
        if SELECT_REPRESENTATIVES_IN_UI:
            non_representative_suids = all_node_suids - representative_suids
            if non_representative_suids:
                print(f"Selecting {len(non_representative_suids)} non-representative nodes in Cytoscape UI...")
                p4c.clear_selection('nodes', network=network_suid)
                p4c.select_nodes(list(non_representative_suids), network=network_suid)
                print("Non-representative nodes selected.")
            else:
                print("No non-representative nodes found to select.")
        
        print("\nScript finished successfully.")
        
    except ImportError as e:
         print(f"Error: Missing required library. {e}")
         print("Please install py4cytoscape: pip install py4cytoscape")
    except Exception as e:
        print(f"\nAn error occurred: {e}")
        print("Please ensure Cytoscape is running, the network is loaded,")
        print("and the specified paths and attributes are correct.")
        import traceback
        traceback.print_exc() # Print detailed traceback for debugging

# --- Run the script ---
if __name__ == '__main__':
    # You can optionally pass a network name or SUID via command line
    # Example: python select_cluster_representatives.py "My Network Name"
    target_network = sys.argv[1] if len(sys.argv) > 1 else None
    select_representatives(network=target_network) 