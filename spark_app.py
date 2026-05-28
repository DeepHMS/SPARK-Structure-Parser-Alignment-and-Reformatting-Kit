import streamlit as st
import pandas as pd
import json
import io
import py3Dmol
from stmol import showmol
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select

# --- Page Configuration ---
st.set_page_config(page_title="SPARK Toolkit", layout="wide", page_icon="⚡")

# --- Helper Classes & Functions ---
class StructureFilter(Select):
    def __init__(self, remove_h=False, remove_hoh=False, allowed_chains=None):
        self.remove_h = remove_h
        self.remove_hoh = remove_hoh
        self.allowed_chains = allowed_chains

    def accept_chain(self, chain):
        if self.allowed_chains is not None and chain.id not in self.allowed_chains:
            return 0
        return 1

    def accept_residue(self, residue):
        if self.remove_hoh and residue.get_resname() in ['HOH', 'WAT']:
            return 0
        return 1

    def accept_atom(self, atom):
        if self.remove_h and atom.element == 'H':
            return 0
        return 1

def get_pdb_stats(structure):
    stats = {'Chains': 0, 'Residues': 0, 'Atoms': 0, 'Hydrogens': 0, 'Water (HOH)': 0}
    chains = []
    for model in structure:
        for chain in model:
            stats['Chains'] += 1
            chains.append(chain.id)
            for residue in chain:
                if residue.get_resname() in ['HOH', 'WAT']:
                    stats['Water (HOH)'] += 1
                else:
                    stats['Residues'] += 1
                for atom in residue:
                    stats['Atoms'] += 1
                    if atom.element == 'H':
                        stats['Hydrogens'] += 1
    return stats, chains

def render_3d_viewer(pdb_string, chain_color, water_color, h_color, bg_color):
    """Renders a single py3Dmol viewer with customized colors and forced Hydrogen rendering."""
    view = py3Dmol.view(width=800, height=500)
    
    view.addModel(pdb_string, 'pdb', {'keepH': True})
    
    # 1. Style Chains
    if chain_color == 'spectrum':
        view.setStyle({'model': -1}, {"cartoon": {'color': 'spectrum'}})
    elif chain_color == 'ssPyMOL':
        view.setStyle({'model': -1}, {"cartoon": {'colorscheme': 'ssPyMOL'}})
    else:
        view.setStyle({'model': -1}, {"cartoon": {'color': chain_color}})
        
    # 2. Style Water and Hydrogens based on user color pickers
    view.addStyle({'resn': 'HOH'}, {'sphere': {'color': water_color, 'radius': 0.5}})
    view.addStyle({'resn': 'WAT'}, {'sphere': {'color': water_color, 'radius': 0.5}})
    view.addStyle({'elem': 'H'}, {'sphere': {'color': h_color, 'radius': 0.4}})
    
    view.setBackgroundColor(bg_color)
    view.zoomTo()
    return view

# --- Main App UI ---
st.title("⚡ SPARK: Structure Parser, Alignment, and Reformatting Kit")
st.markdown("SPARK is a lightweight Python toolkit bridging classic PDB with modern AI models. Seamlessly convert mmCIF & JSON formats to PDB, parse complex coordinate data, and run high-precision structural alignments to compare and evaluate multiple models of the same protein.")

# ==========================================
# GLOBAL SIDEBAR CONTROLS
# ==========================================
with st.sidebar:
    st.header("🎨 Advanced Color Settings")
    st.markdown("These colors will update all 3D viewers across every tab.")
    
    chain_color = st.selectbox("Chain Color Scheme", ['spectrum', 'ssPyMOL', 'blue', 'green', 'red', 'orange', 'purple', 'gray', 'white'])
    water_color = st.color_picker("Water (HOH) Color", "#FF0000") # Default Red
    h_color = st.color_picker("Hydrogen (H) Color", "#D3D3D3")     # Default Light Gray
    bg_color = st.color_picker("Background Color", "#FFFFFF")     # Default White
    
    st.markdown("---")
    
    st.header("🎮 Camera & Export Controls")
    st.markdown("""
    **Interact directly with the 3D Viewer:**
    - 🔍 **Zoom In/Out:** Use Mouse Scroll Wheel
    - 🔄 **Rotate:** Left Click + Drag
    - ↔️ **Pan Axis:** Right Click + Drag *(or Ctrl + Left Click + Drag)*
    - 📸 **Screenshot:** Right Click on the 3D model ➔ Select *'Save image as...'*. This captures a high-res PNG exactly as you see it!
    """)
    
    if st.button("🔄 Reset Viewer Angles"):
        st.rerun()

# ==========================================
# Create Tabs
# ==========================================
tab1, tab2, tab3 = st.tabs([
    "1️⃣ Format Conversions", 
    "2️⃣ Preprocessing & Clean-Up", 
    "3️⃣ Structure Model Comparison"
])

# ==========================================
# TAB 1: Format Conversions
# ==========================================
with tab1:
    st.header("Format Conversions (JSON / mmCIF ➡️ PDB)")
    uploaded_file = st.file_uploader("Upload .json or .cif file", type=['json', 'cif'], key="t1_uploader")
    
    if uploaded_file is not None:
        file_ext = uploaded_file.name.split('.')[-1].lower()
        extracted_pdb_string = None
        
        try:
            if file_ext == 'json':
                data = json.load(uploaded_file)
                if "structures" in data:
                    raw_text = data["structures"][0]["structure"]
                elif "structures_in_ranked_order" in data:
                    raw_text = data["structures_in_ranked_order"][0]["structure"]
                else:
                    raw_text = None
                
                if raw_text:
                    if "_audit_conform.dict_name" in raw_text:
                        cif_io = io.StringIO(raw_text)
                        parser = MMCIFParser(QUIET=True)
                        structure = parser.get_structure("temp", cif_io)
                        pdb_io = PDBIO()
                        pdb_io.set_structure(structure)
                        out_io = io.StringIO()
                        pdb_io.save(out_io)
                        extracted_pdb_string = out_io.getvalue()
                    else:
                        extracted_pdb_string = raw_text
            
            elif file_ext == 'cif':
                cif_text = uploaded_file.getvalue().decode("utf-8")
                cif_io = io.StringIO(cif_text)
                parser = MMCIFParser(QUIET=True)
                structure = parser.get_structure("temp", cif_io)
                pdb_io = PDBIO()
                pdb_io.set_structure(structure)
                out_io = io.StringIO()
                pdb_io.save(out_io)
                extracted_pdb_string = out_io.getvalue()

            if extracted_pdb_string:
                col1, col2 = st.columns([1, 2])
                with col1:
                    st.download_button("⬇️ Download Converted PDB", data=extracted_pdb_string, file_name=f"{uploaded_file.name.split('.')[0]}_converted.pdb", mime="chemical/x-pdb")
                with col2:
                    view = render_3d_viewer(extracted_pdb_string, chain_color, water_color, h_color, bg_color)
                    showmol(view, height=500, width=800)
                    
        except Exception as e:
            st.error(f"Error during conversion: {e}")

# ==========================================
# TAB 2: Structure Preprocessing & Clean-Up
# ==========================================
with tab2:
    st.header("Structure Preprocessing & Clean-Up")
    raw_pdb = st.file_uploader("Upload a Single .pdb Structure", type=['pdb'], key="t2_uploader")
    
    if raw_pdb:
        pdb_string = raw_pdb.getvalue().decode("utf-8")
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", io.StringIO(pdb_string))
        
        init_stats, available_chains = get_pdb_stats(structure)
        
        col_ctrl, col_view = st.columns([1, 2])
        
        with col_ctrl:
            st.subheader("Filter Settings")
            remove_h = st.checkbox("Remove Hydrogens (H)", value=False)
            remove_hoh = st.checkbox("Remove Water (HOH/WAT)", value=False)
            selected_chains = st.multiselect("Select Chains to Keep", options=list(set(available_chains)), default=list(set(available_chains)))
            
            io_writer = PDBIO()
            io_writer.set_structure(structure)
            filtered_io = io.StringIO()
            io_writer.save(filtered_io, select=StructureFilter(remove_h=remove_h, remove_hoh=remove_hoh, allowed_chains=selected_chains))
            
            filtered_pdb_string = filtered_io.getvalue()
            filtered_struct = parser.get_structure("filtered", io.StringIO(filtered_pdb_string))
            final_stats, _ = get_pdb_stats(filtered_struct)
            
            st.subheader("Statistics")
            stat_df = pd.DataFrame([init_stats, final_stats], index=["Original", "Filtered"]).T
            st.dataframe(stat_df, use_container_width=True)
            
            st.download_button("⬇️ Download Filtered PDB", data=filtered_pdb_string, file_name=f"filtered_{raw_pdb.name}", mime="chemical/x-pdb")
            
        with col_view:
            st.subheader("Interactive Viewer")
            view2 = render_3d_viewer(filtered_pdb_string, chain_color, water_color, h_color, bg_color)
            showmol(view2, height=600, width=800)

# ==========================================
# TAB 3: Structure Model Comparison
# ==========================================
with tab3:
    st.header("Structure Model Comparison & Synchronized Viewer")
    num_files = st.number_input("How many PDB structures do you want to compare?", min_value=2, max_value=12, value=3)
    
    pdb_data_list = []
    labels = []
    stats_list = []
    
    upload_cols = st.columns(3)
    for i in range(num_files):
        col_idx = i % 3
        with upload_cols[col_idx]:
            label = st.text_input(f"Label for Model {i+1}", value=f"Model {i+1}", key=f"label_{i}")
            file = st.file_uploader(f"Upload PDB {i+1}", type=['pdb'], key=f"file_{i}")
            
            if file:
                pdb_str = file.getvalue().decode("utf-8")
                pdb_data_list.append(pdb_str)
                labels.append(label)
                
                parser = PDBParser(QUIET=True)
                struct = parser.get_structure(label, io.StringIO(pdb_str))
                stats, _ = get_pdb_stats(struct)
                stats['Model Name'] = label
                stats_list.append(stats)
    
    if len(pdb_data_list) > 0:
        st.markdown("---")
        st.subheader("📊 Comparative Statistics")
        df_compare = pd.DataFrame(stats_list).set_index('Model Name')
        st.dataframe(df_compare, use_container_width=True)
        
        st.subheader("🔬 Interactive Grid Viewer")
        
        # --- NEW SYNCHRONIZATION CHECKBOX ---
        sync_views = st.checkbox("🔄 Synchronize Viewers (Zoom/Pan)", value=True, help="If checked, interacting with one model will move all models. Uncheck to move them independently.")
        
        rows = (len(pdb_data_list) + 2) // 3
        
        # Apply the sync_views boolean to the linked parameter
        grid_view = py3Dmol.view(viewergrid=(rows, 3), width=1200, height=400 * rows, linked=sync_views)
        grid_view.setBackgroundColor(bg_color)
        
        for i, pdb_str in enumerate(pdb_data_list):
            r, c = divmod(i, 3)
            
            grid_view.addModel(pdb_str, 'pdb', {'keepH': True}, viewer=(r, c))
            
            if chain_color == 'spectrum':
                grid_view.setStyle({'model': -1}, {"cartoon": {'color': 'spectrum'}}, viewer=(r, c))
            elif chain_color == 'ssPyMOL':
                grid_view.setStyle({'model': -1}, {"cartoon": {'colorscheme': 'ssPyMOL'}}, viewer=(r, c))
            else:
                grid_view.setStyle({'model': -1}, {"cartoon": {'color': chain_color}}, viewer=(r, c))
                
            grid_view.addStyle({'resn': 'HOH'}, {'sphere': {'color': water_color, 'radius': 0.5}}, viewer=(r, c))
            grid_view.addStyle({'resn': 'WAT'}, {'sphere': {'color': water_color, 'radius': 0.5}}, viewer=(r, c))
            grid_view.addStyle({'elem': 'H'}, {'sphere': {'color': h_color, 'radius': 0.4}}, viewer=(r, c))
            
            # Label contrast based on background
            font_col = 'black' if bg_color.upper() == '#FFFFFF' else 'white'
            grid_view.addLabel(labels[i], {'position': {'x':0, 'y':0, 'z':0}, 'useScreen': True, 'fontColor': font_col, 'backgroundColor': 'transparent', 'fontSize': 16, 'alignment': 'topLeft'}, viewer=(r, c))
            grid_view.zoomTo(viewer=(r, c))
            
        showmol(grid_view, height=400 * rows, width=1200)
