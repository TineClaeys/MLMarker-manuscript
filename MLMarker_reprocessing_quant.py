import logging
import glob
import pandas as pd
import numpy as np
import os
import re
import warnings
# Filter out the specific warning message
warnings.filterwarnings("ignore", message="A value is trying to be set on a copy of a DataFrame or Series through chained assignment using an inplace method.")

import plotly.graph_objects as go
import sys

# import sys
# sys.modules.pop("mlmarker", None)
import mlmarker
# import mlmarker
# print(mlmarker.__file__)
import mlmarker.model as MLMarker

logging.FileHandler('MLMarker.log') #Create a log file
logging.basicConfig(level=logging.DEBUG, handlers=[logging.FileHandler('MLMarker.log')])


def extract_pxd(file_path):
    match = re.search(r'PXD.*?/(?<!zip/)', file_path) # Pattern that starts with PXD and ends with / but not with zip/
    if match:
        return match.group(0)[:-1]
    else:
        return None


def extract_filename(file_path):
    match = re.search(r'/(?!PXD).*?/', file_path) # Pattern that starts with / and everything else other than PXD and ends with /
    if match:
        return match.group(0)[1:-1]
    else:
        return None


    
def ionbot_first_to_prot_bin(df, num_psms=False):
        """Only search in the rows : q-value <= 0.01, database = T, proteins != CONTAMINANT and '|'  
        Return a list of uniprot id"""
        df = df[['proteins', 'database', 'q-value']]
        df = df[(df['database']=='T') & (df['q-value']<=0.01)]
        if num_psms:
            num_psms = df.shape[0]
        df = df[~df['proteins'].apply(lambda x: '|' in x or 'CONTAMINANT' in x)]
        uniprot_df = df['proteins'].apply(lambda uniprot: re.search(r'\)\)\(\((.*?)\)\)', uniprot).group(1) 
                                          if re.search(r'\)\)\(\((.*?)\)\)', uniprot) else None) # Found the uniprot id with a ))((.*)) pattern
        uniprot_df = uniprot_df.dropna() #Maybe not necessary 
        list_uniprot = (list(set(uniprot_df)))
        if num_psms:
            return list_uniprot, num_psms
        else:
            return list_uniprot

def ionbot_first_to_prot_quant(df, lengths):
        """Only search in the rows : q-value <= 0.01, database = T, proteins != CONTAMINANT and '|'  
        Return a list of uniprot id"""
        logging.info("Starting ionbot_first_to_prot_quant")
        if 'proteins' not in df.columns:
            logging.error("Column 'proteins' is missing in the input DataFrame. Terminating function.")
            return pd.DataFrame()
        df = df[['proteins', 'database', 'q-value',]]
        df = df[(df['database']=='T') & (df['q-value']<=0.01)]
        df = df[~df['proteins'].apply(lambda x: '|' in x or 'CONTAMINANT' in x)]
        #if nothing remains after filtering, return an empty DataFrame
        if df.shape[0] == 0:
            logging.error("No valid proteins found after filtering. Terminating function.")
            return pd.DataFrame()
        logging.info(df)
        uniprot_df = df['proteins'].apply(lambda uniprot: re.search(r'\)\)\(\((.*?)\)\)', uniprot).group(1) 
                                          if re.search(r'\)\)\(\((.*?)\)\)', uniprot) else None) # Found the uniprot id with a ))((.*)) pattern
        uniprot_df = uniprot_df.value_counts().to_frame().reset_index()
        uniprot_df = uniprot_df.merge(lengths, left_on='proteins', right_on='Entry')
        uniprot_df = uniprot_df[uniprot_df['Length'] != np.nan]
        uniprot_df['count'] = uniprot_df['count'].astype(float)
        uniprot_df['Length'] = uniprot_df['Length'].astype(float)
        uniprot_df['SAF'] = uniprot_df['count']/uniprot_df['Length']
        total_SAF = uniprot_df['SAF'].sum()
        uniprot_df['NSAF'] = uniprot_df['SAF']/total_SAF
        logging.info("Done with quant")
        return uniprot_df

def ionbot_first_to_mods(df):
    """Only search in the rows : q-value <= 0.01, database = T, proteins != CONTAMINANT and '|'  
    Return a list of uniprot id"""
    df = df[['proteins', 'database', 'q-value', 'modifications']]
    df = df[(df['database']=='T') & (df['q-value']<=0.01)]
    df = df[~df['proteins'].apply(lambda x: '|' in x or 'CONTAMINANT' in x)]
    uniprot_df = df['modifications'].str.split('_or_').explode()
    uniprot_df = uniprot_df.value_counts().to_frame().reset_index()
    return uniprot_df

def csv_path(path_file):   
    """Return a dictionnary with the path of the csv files"""
    for result_dir in glob.iglob(path_file, recursive=True): #Recursive search of the zip files in the folder
        pxd = os.path.basename(os.path.split(result_dir)[0])
        for run_dir in os.listdir(result_dir):
            if run_dir.endswith('.mgf.gzip'): #Check if the file is a mgf.gz file
                result_file = os.path.join(result_dir, run_dir, "ionbot.first.csv")
                if os.path.exists(result_file):
                    yield pxd, run_dir, result_file


def MLMarker_reprocessing_binary(prot_list):
    data = np.ones((1, len(prot_list)))
    predict_df = pd.DataFrame(data, columns=prot_list)
    model = mlmarker.MLMarker(binary=True, dev=dev, penalty_factor=1)
    model.load_sample(predict_df.iloc[0:1,:])
    prediction = model.explainability.adjusted_shap_values_df(penalty_factor=1, n_preds=100)
    prediction = prediction.sum(axis=1).sort_values(ascending=False)
    df_predi = pd.DataFrame(prediction).transpose()
    return df_predi

def csv_path(pxd, base_path):
    """Generate paths to the csv files for a given PXD."""
    search_path = os.path.join(base_path, f"{pxd}/IONBOT_v0.11.3")
    for result_dir in glob.iglob(search_path, recursive=True):
        for run_dir in os.listdir(result_dir):
            if run_dir.endswith('.mgf.gzip'):
                result_file = os.path.join(result_dir, run_dir, "ionbot.first.csv")
                if os.path.exists(result_file):
                    yield pxd, run_dir, result_file

def MLMarker_to_Dataframe(pxd, base_path, lengths, dev, tissue_level, count_psms=False):
    """Return a dataframe with the data of the mgf.gz files after the MLMarker analysis"""
    logging.info("Starting MLMarker_to_Dataframe")
    df_big = pd.DataFrame()
    processed = set()
    prot_dict, psm_dict = {}, {}
    for pxd, mgf_file, result_file in csv_path(pxd, base_path):
        logging.info(f"Now running {result_file}")
        if (pxd, mgf_file) in processed:
            logging.error(f'Duplicate in {pxd}/{mgf_file}')
            continue
        
        if not result_file.endswith('ionbot.first.csv'):
            continue

        try:

            with open(result_file) as f: 
                df_first = pd.read_csv(f)
                #log filename
                logging.info(f"started processing {mgf_file}")
                logging.info(df_first.shape)
                 # Early exit if 'proteins' column is not in the DataFrame
                if 'proteins' not in df_first.columns:
                    logging.error(f"'proteins' column missing in {pxd}/{mgf_file}. Skipping this file.")
                    continue  # Skip this file and move to the next

                if df_first.shape[0] == 0:
                    logging.error(f'Empty file in {pxd}/{mgf_file}')
                    continue
        
        except Exception as e:
            # Handle file reading or parsing errors
            logging.exception(f"Error reading or processing the result file for {pxd}/{mgf_file}. Skipping to the next file.")
            continue
            
        prot_list = ionbot_first_to_prot_quant(df_first, lengths)
        if prot_list.empty:
            logging.error(f'No valid proteins found in {pxd}/{mgf_file}. Skipping to the next file.')
            continue
        if count_psms:
            unique_proteins, num_psms = ionbot_first_to_prot_bin(df_first, num_psms=count_psms)
            psm_dict[mgf_file] = num_psms
        else:
            unique_proteins = ionbot_first_to_prot_bin(df_first)
        
        prot_dict[mgf_file] = unique_proteins

        logging.info(f'starting prediction of {mgf_file}')
        df_temp = MLMarker_reprocessing_quant(prot_list, dev, tissue_level)
        logging.info('Done with predicting')
        logging.info(df_temp.head())
        df_temp['PXD_folder'] = pxd
        df_temp['filename'] = mgf_file

        df_big = pd.concat([df_big, df_temp])
        processed.add((pxd, mgf_file))

    if count_psms:
        return df_big, prot_dict, psm_dict
    else:
        return df_big, prot_dict



def MLMarker_reprocessing_quant(prot_list, dev, tissue_level):
    logging.info("Starting MLMarker_reprocessing_quant")
    
    data = prot_list.pivot_table(columns='proteins', values='NSAF', aggfunc='sum')
    data = data.fillna(0)
    if data.empty:
        raise ValueError("Input data is empty after pivoting. Check `prot_list`.")
    mlmarker_model = mlmarker.MLMarker(binary=False, dev=dev, penalty_factor = 1)
    mlmarker_model.load_sample(data.iloc[0:1, :])
    added_features = mlmarker_model.load_sample(data.iloc[0:1, :], output_added_features=True)
    logging.info(len(added_features))
    # Get model features
    features = mlmarker_model.get_model_features()
    logging.info(len(features))
    # Compute missing features penalty
    missing_features = set(features) - set(data.columns.tolist())
    #penalty_factor = len(missing_features) / len(features) * 100
    # penalty_factor = 1
    # # Create a new MLMarker instance with penalty factor
    # mlmarker_model.penalty_factor = penalty_factor
    
    prediction = mlmarker_model.explainability.adjusted_absent_shap_values_df(n_preds=100)

    if tissue_level:
        prediction = prediction.sum(axis=1).sort_values(ascending=False)
    df_predi = pd.DataFrame(prediction).transpose()
    return df_predi


# Main function to call in the notebook
def get_mlmarker_dataframe(pxd, base_path, lengths_path, dev, tissue_level=True, count_psms=False):
    """
    Given a PXD identifier, return the MLMarker dataframe.
    
    Args:
        pxd (str): The PXD identifier.
        base_path (str): Base directory where PRIDE data is stored.
        lengths_path (str): Path to the uniprot lengths file.

    Returns:
        pd.DataFrame: Processed MLMarker dataframe.
    """
    lengths = pd.read_csv(lengths_path, sep='\t')
    model = mlmarker.MLMarker(binary=False, dev=dev, penalty_factor=0)
    if count_psms:
        mlmarker_df, unique_protein_dict, psm_dict = MLMarker_to_Dataframe(pxd, base_path, lengths, dev, tissue_level, count_psms)
        return mlmarker_df, unique_protein_dict, psm_dict
    else:
        mlmarker_df, unique_protein_dict = MLMarker_to_Dataframe(pxd, base_path, lengths, dev, tissue_level, count_psms)
        return mlmarker_df, unique_protein_dict



def get_protein_info(protein_id):
    """
    Get protein information from UniProt.
    
    Parameters:
        protein_id (str): UniProt protein ID.
    
    Returns:
        dict: Protein information.
    """
    import bioservices
    u = bioservices.UniProt()
    try:
        protein_info = u.search(protein_id, columns="accession, id, protein_name, cc_tissue_specificity")
        protein_info = protein_info.split('\n')[1].split('\t')
        protein_dict = {
            'id': protein_info[0],
            'entry name': protein_info[1],
            'protein_names': protein_info[2]
        }
        if len(protein_info) == 4:
            protein_dict['tissue_specificity'] = protein_info[3]
        return protein_dict
    except:
        print(f"Error retrieving information for protein {protein_id}")
        return None

def get_go_enrichment(protein_list):
    from gprofiler import GProfiler
    import plotly.graph_objects as go
    # Initialize g:Profiler
    gp = GProfiler(return_dataframe=True)

    # Dictionary to store GO terms and p-values for each tissue
    go_dict = {}


    # Perform GO enrichment
    results = gp.profile(organism='hsapiens', query=protein_list, sources=['GO:BP', 'GO:MF', 'GO:CC', 'HPA'], combined=True)
    results = results[results['p_value']< 0.05]
    # Store results in the dictionary: {tissue: {GO_term: p-value}}
    return results

def visualise_custom_plot(df):
        
    # Aggregate positive and negative contributions per tissue
    positive_totals = df.clip(lower=0).sum(axis=1)
    negative_totals = df.clip(upper=0).abs().sum(axis=1)

    # Create the figure
    fig = go.Figure()

    # Add positive contributions (green bars)
    fig.add_trace(
        go.Bar(
            x=df.index,
            y=positive_totals,
            name="Positive Contributions",
            marker_color='green',
            hoverinfo='x+y',
        )
    )

    # Add negative contributions (red bars)
    fig.add_trace(
        go.Bar(
            x=df.index,
            y=negative_totals,
            name="Negative Contributions",
            marker_color='red',
            hoverinfo='x+y',
        )
    )

    # Customizing layout
    fig.update_layout(
        barmode='group',  # Group positive and negative bars side-by-side
        title='Grouped Barplot of Total Protein Contributions by Tissue',
        xaxis_title='Tissues',
        yaxis_title='Total Contributions',
        xaxis=dict(tickangle=-45),  # Tilt the x-axis labels for better readability
        template="plotly_white"
    )

    fig.show()

import plotly.graph_objects as go

def visualise_custom_tissue_plot(df, tissue_name, top_n=10, show_others=False, threshold_others = 0.001):
    df = df.loc[[tissue_name]]

    # Separate positive and negative values for the tissue
    positive_contributions = df.clip(lower=0)  # Keep only positive values
    negative_contributions = df.clip(upper=0).abs()  # Keep absolute values of negatives

    # Filter significant contributions
    positive_main = positive_contributions.loc[:, (positive_contributions > threshold_others).any()]
    positive_others = positive_contributions.loc[:, (positive_contributions <= threshold_others).all()].sum(axis=1)

    negative_main = negative_contributions.loc[:, (negative_contributions > threshold_others).any()]
    negative_others = negative_contributions.loc[:, (negative_contributions <= threshold_others).all()].sum(axis=1)

    # Sort positive and negative contributions by total value
    sorted_positive = positive_main.sum(axis=0).sort_values(ascending=False)
    sorted_negative = negative_main.sum(axis=0).sort_values(ascending=False)

    # Select top N positive and negative proteins
    top_positive_contributions = sorted_positive.head(top_n).index.tolist()
    top_negative_contributions = sorted_negative.head(top_n).index.tolist()

    # Plotting
    fig = go.Figure()

    # Add all positive contributions (green bars)
    for protein in sorted_positive.index:
        # Check if the protein is one of the top N and add its label
        is_top = protein in top_positive_contributions
        fig.add_trace(
            go.Bar(
                x=positive_contributions.index,
                y=positive_main[protein],
                name=protein,
                marker_color="green" if is_top else "darkgreen",
                hoverinfo="name+y",
                hoverlabel=dict(namelength=-1),
                showlegend=False,
                text=protein if is_top else None,  # Show label for top proteins
                textposition="outside",  # Position the label inside the bar
                cliponaxis=False,  # Allow the label to be outside the bar
            )
        )
    # Add lines for top proteins to connect labels outside the bars
    for protein in top_positive_contributions:
        fig.add_trace(
            go.Scatter(
                x=[positive_contributions.index[0], positive_contributions.index[0]],
                y=[positive_contributions[protein].min(), positive_contributions[protein].max()],
                mode="lines+text",
                line=dict(color="green", width=2, dash="dot"),  # Line connecting label to bar
                text=[protein],
                textposition="middle right",
                showlegend=False,
                textfont=dict(color="green", size=12)
            )
        )
    # Add "Others" for positive contributions
    if show_others and positive_others.sum() > 0:
        fig.add_trace(
            go.Bar(
                x=positive_contributions.index,
                y=positive_others,
                name="Others (Positive)",
                marker_color="lightgreen",
                hoverinfo="name+y",
                hoverlabel=dict(namelength=-1),
                showlegend=False,
            )
        )

  # Add negative contributions (sorted by total contribution)
    for protein in sorted_negative.index:
        is_top = protein in top_negative_contributions
        fig.add_trace(
            go.Bar(
                x=negative_contributions.index,
                y=negative_main[protein],
                name=protein,
                marker_color="red" if is_top else "darkred",
                hoverinfo="name+y",
                hoverlabel=dict(namelength=-1),
                showlegend=False,
                text=protein if is_top else None,  # Show label for top proteins
                textposition="outside",  # Position the label outside the bar
                cliponaxis=False,  # Allow the label to be outside the bar
            )
        )

    # Add "Others" for negative contributions
    if show_others and negative_others.sum() > 0:
        fig.add_trace(
            go.Bar(
                x=negative_contributions.index,
                y=negative_others,
                name="Others (Negative)",
                marker_color="lightcoral",
                hoverinfo="name+y",
                hoverlabel=dict(namelength=-1),
                showlegend=False,
            )
        )

    # Customizing layout
    fig.update_layout(
        barmode="stack",  # Stack the bars
        title=f"""Protein Contributions for {tissue_name} (threshold={threshold_others})""",
        xaxis_title="Cluster",
        yaxis_title="Protein Contributions",
        xaxis={"categoryorder": "array", "categoryarray": sorted_positive.index.tolist() + sorted_negative.index.tolist()},
        hovermode="closest",
        template="plotly_white",
        width=600,
        height=800,
        margin=dict(l=100, r=100),  # Adjust margins
    )
    fig.show()

def prediction_df_2tissues_scatterplot(df, tissues=list):
    tissueA = tissues[0]
    tissueB = tissues[1]
    df_vis = df.T
    fig = go.Figure(data=go.Scatter(
        x=df_vis[tissueA],
        y=df_vis[tissueB],
        mode='markers',
        marker=dict(
            size=8,
            color='blue',  # You can change the color here
            opacity=0.7
        ),
        text=[f"Protein: {protein}<br>{tissueA} SHAP: {pg_shap}<br>{tissueB} value: {brain_value}" 
            for protein, pg_shap, brain_value in zip(df_vis.index, df_vis[tissueA], df_vis[tissueB])],
        hoverinfo='text'
    ))

    fig.update_layout(
        title=f'Scatterplot of {tissueA} SHAP values vs {tissueB} values',
        xaxis_title=f'{tissueA} SHAP values',
        yaxis_title=f'{tissueB} SHAP values',
        xaxis=dict(color='black', zeroline=True, zerolinecolor='darkgrey'),
        yaxis=dict(color='black', zeroline=True, zerolinecolor='darkgrey')
    )

    fig.show()