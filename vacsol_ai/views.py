import pickle
import shutil
import subprocess
import zipfile
import iedb
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from PyBioMed.PyPretreat import PyPretreatPro
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from PyBioMed import Pyprotein
import pandas as pd
import os
import fastaparser
from sklearn import preprocessing as per
import sys    
from django.http import HttpRequest, HttpResponse, HttpResponseRedirect, JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.shortcuts import redirect, render
import pickle
import pandas as pd
import os
import sys
import logging
import biolib
from sklearn.feature_selection import VarianceThreshold

@csrf_exempt
def home(request):
    return render (request, 'index.html')

def get_latest_logs(request):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    log_file_path = os.path.join(script_dir, "vacsol_ai.log")
    with open(log_file_path, 'r') as log_file:
        latest_logs = log_file.readlines()
    return JsonResponse({'logs': latest_logs})

def progress_callback(message, progress):
    logger = logging.getLogger(__name__)
    log_message = f"{message}"
    logger.info(log_message, extra={'progress': progress})

### ESKAPE ###

def calculate_features(file_path, progress_callback):
    total_steps = 8
    completed_steps = 0
    # Use os.path.abspath to get the absolute file path
    script_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(script_dir, "sequences.fasta")
    file_path = os.path.abspath(file_path)
    MERGED_DATAFRAMES = []
    with open(file_path) as fasta_file:
        parser = fastaparser.Reader(fasta_file)
        seq_num = 1
        ID_list = []
        list_data = []
        for seq in parser:
            ID = seq.id
            sequence = seq.formatted_fasta()
            sequence1 = seq.sequence_as_string()
            ID_list.append(ID)
            df_ID = pd.DataFrame(ID_list, columns=['ID'])
            completed_steps += 1
            progress=int((completed_steps / total_steps) * 100)
            progress_callback(f"Step 1: Sequence {ID} uploaded successfully", 1)

            # Send POST request to MHC class I peptide binding prediction tool:
            mhci_res1 = iedb.query_mhci_binding(method="recommended", sequence= sequence, allele="HLA-A*01:01", length="9")
            mhci_res2 = iedb.query_mhci_binding(method="recommended", sequence= sequence, allele="HLA-A*02:01", length="9")
            #concat mhci-res
            mhci_res = pd.concat([mhci_res1, mhci_res2], axis=0)
            
            # Send POST request to MHC class II peptide binding prediction tool:
            mhcii_res1 = iedb.query_mhcii_binding(method="recommended", sequence= sequence, allele="HLA-DRB1*01:01", length=None)
            mhcii_res2 = iedb.query_mhcii_binding(method="recommended", sequence= sequence, allele="HLA-DRB1*01:02", length=None)
            
            #concat mhcii-res
            mhcii_res = pd.concat([mhcii_res1, mhcii_res2], axis=0)
            
            # Send POST request to B-cell epitope prediction tool:
            bcell_res = iedb.query_bcell_epitope(method="Bepipred", sequence= sequence1, window_size=9)
            completed_steps += 1
            progress=int((completed_steps / total_steps) * 100)
            progress_callback("Step 2: Epitope Analysis Completed", 2)

            # Send POST request to surface probability prediction tool:
            sprob_res = iedb.query_bcell_epitope(method="Emini", sequence= sequence1, window_size=9)
            completed_steps += 1
            progress=int((completed_steps / total_steps) * 100)
            progress_callback("Step 3: Adhesion Probability Analysis Completed", 3)

            # Send POST request to antigenicity prediction tool:
            antigenicity_1 = iedb.query_bcell_epitope(method="Kolaskar-Tongaonkar", sequence= sequence1, window_size=9)
            completed_steps += 1
            progress=int((completed_steps / total_steps) * 100)
            progress_callback("Step 4: Antigenicity Prediction Completed", 4)

            # Getting means - mhci
            mhci_res['score'] = pd.to_numeric(mhci_res['score'], errors='coerce')
            mhci_res['percentile_rank'] = pd.to_numeric(mhci_res['percentile_rank'], errors='coerce')

            df_score_mhci = mhci_res["score"].mean()
            df_rank_mhci = mhci_res["percentile_rank"].mean()

            # Getting means - mhcii
            mhcii_res['score'] = pd.to_numeric(mhcii_res['score'], errors='coerce')
            mhcii_res['rank'] = pd.to_numeric(mhcii_res['rank'], errors='coerce')

            rank_mhcii = mhcii_res["rank"].mean()
            score_mhcii = mhcii_res["score"].mean()

            # Getting means - bcells
            bcell_res['Score'] = pd.to_numeric(bcell_res['Score'], errors='coerce')
            df_bcells_final = bcell_res["Score"].mean()

            # Getting means - s_probability
            sprob_res['Score'] = pd.to_numeric(sprob_res['Score'], errors='coerce')
            df_sprob_final = sprob_res["Score"].mean()

            # Getting means - antigenicity(scale1)
            antigenicity_1['Score'] = pd.to_numeric(antigenicity_1['Score'], errors='coerce')
            antigenicity_1_final = antigenicity_1["Score"].mean()

            #features_file for users
            script_dir = sys.path[0]
            results_folder_path = os.path.join(script_dir, "Analysis_Results")

            #Epitope_analysis and Antigenicity score
            epitope_filename = 'epitope_analysis.csv'
            column_names_ii = ['Protein_ID', 'antigenicity', 'b_cells_probability_score', 'mhci_probability_score', 'mhci_rank', 'mhcii_rank', 'mhcii_score']
            epitope_data = [ID, antigenicity_1_final, df_bcells_final, df_score_mhci, df_rank_mhci, rank_mhcii, score_mhcii]
            df_epitope = pd.DataFrame(data=[epitope_data], columns=column_names_ii)
            epitope_result_path = os.path.join(results_folder_path, epitope_filename)
            df_epitope.to_csv(epitope_result_path, index=False)

            #adhesion_probability score
            adhesion_filename = 'adhesion_probability.csv'
            column_names_iii = ['Protein_ID' , 'adhesion_probability']
            adhesion_data = [ID, df_sprob_final]
            df = pd.DataFrame(data=[adhesion_data], columns=column_names_iii)
            df_adhesion = pd.DataFrame(data=[adhesion_data], columns=column_names_iii)
            adhesion_result_path = os.path.join(results_folder_path, adhesion_filename)
            df_adhesion.to_csv(adhesion_result_path, index=False)
            
            #Analysis of Physiochemical Features:
            X = ProteinAnalysis(sequence1)
            #by protoparam
            scale1= "%0.2f" % X.molecular_weight()
            scale2= "%0.2f" % X.charge_at_pH(7)
            scale3= "%0.2f" % X.isoelectric_point()
            scale4= "%0.2f" % X.gravy()
            scale5= "%0.2f" % X.instability_index()
            scale6= "%0.2f" % X.aromaticity()
            #by pybiomed
            aa_no= PyPretreatPro.ProteinCheck(sequence1)
            
            #Final CSV
            header = ['Protein_ID','aa_no','molecular_weight', 'pH_charge', 'isoelectric point', 'Hydropathy_gravy', 'instability_index', 'aromaticity', 'antigenicity_1', 'b_cells_probability_score', 'mhci_probability_score', 'mhci_rank', 'mhcii_rank', 'mhcii_score', 'surface_probability']
            data = [ID, aa_no, scale1, scale2, scale3, scale4, scale5, scale6, antigenicity_1_final, df_bcells_final, df_score_mhci, df_rank_mhci, rank_mhcii, score_mhcii, df_sprob_final]
            df = pd.DataFrame(data=[data], columns=header)
            
            #print(df)
            list_data.append(df)

            #by pybiomed
            aa_no= PyPretreatPro.ProteinCheck(sequence1)
            protein_class = Pyprotein.PyProtein(sequence1)
            descriptor1 = protein_class.GetCTD()
            descriptor2 = protein_class.GetGearyAuto()
            descriptor3 = protein_class.GetMoranAuto()
            descriptor4 = protein_class.GetMoreauBrotoAuto()
            completed_steps += 1
            progress=int((completed_steps / total_steps) * 100)
            progress_callback("Step 5: Physicochemical Analysis completed", 5)
            
            #df1 = pd.DataFrame.from_dict(descriptor1)
            df1 = pd.DataFrame.from_dict(descriptor1,orient='index').transpose()
            df2 = pd.DataFrame.from_dict(descriptor2,orient='index').transpose()
            df3 = pd.DataFrame.from_dict(descriptor3,orient='index').transpose()
            df4 = pd.DataFrame.from_dict(descriptor4,orient='index').transpose()
            
            df_all = [df1, df2, df3, df4]
            con1 = pd.concat(df_all, axis="columns")
            MERGED_DATAFRAMES.append(con1)

    con = pd.concat(list_data).reset_index()
    print(con)

    # concatenate the script directory with the filename to create the full path of the CSV file
    script_dir = sys.path[0]
    final_csv_path = os.path.join(script_dir, "FINAL.csv")
    completed_steps += 1

    con.to_csv(final_csv_path, index=False)

    #prediction of signal peptides (SignalP. version 6.0)

    csv_folder_path = os.path.join(script_dir, "biolib_results")

    signalp_6 = biolib.load('DTU/SignalP_6')
    signalp6_job = signalp_6.cli(args='--fastafile vacsol_ai/sequences.fasta --output_dir output')
    signalp6_job.save_files('biolib_results')

    #subprocess.run('biolib run DTU/SignalP_6 --fastafile vacsol_ai/sequences.fasta --output_dir output')

    #read signalp results:    
    sp_table_path = os.path.join(csv_folder_path, "prediction_results.txt")
    df = pd.read_table(sp_table_path, sep="\t", header=None, skiprows=1)
    df.columns = ["ID", "Prediction", "OTHER", "SP(Sec/SPI)", "LIPO(Sec/SPII)", "TAT(Tat/SPI)", "TATLIPO(Sec/SPII)", "PILIN(Sec/SPIII)", "CS Position"]
    df.drop(index=df.index[0], axis=1, inplace=True)

    #df.to_csv('E:/ASAB/VacSol-AI/VacSol-ML-ESKAPE-/signal_peptide_analysis.csv', index=False)
    
    #extract the scores
    SP = df['SP(Sec/SPI)'].values
    LIPO = df['LIPO(Sec/SPII)'].values
    TAT = df['TAT(Tat/SPI)'].values
    TATLIPO = df['TATLIPO(Sec/SPII)'].values
    PILIN = df['PILIN(Sec/SPIII)'].values
    OTHER = df['OTHER'].values
    
    completed_steps += 1
    progress=int((completed_steps / total_steps) * 100)
    progress_callback("Step 6: Signal Peptide Analysis Completed", 6)
    add_scores = pd.read_csv(final_csv_path)
    add_scores["signal_peptide_SP"] = SP
    add_scores["signal_peptide_LIPO"] = LIPO
    add_scores["signal_peptide_TAT"] = TAT
    add_scores["signal_peptide_TATLIPO"] = TATLIPO
    add_scores["signal_peptide_PILIN"] = PILIN
    add_scores["signal_peptide_OTHER"] = OTHER

    MERGED = pd.concat(MERGED_DATAFRAMES).reset_index()
    completed_steps += 1
    progress=int((completed_steps / total_steps) * 100)
    progress_callback("Step 7: Creating Dataset", 7)

    data_CSV = [add_scores, MERGED]
    final_con = pd.concat(data_CSV, axis="columns").reset_index()
    #FINAL_CSV.append(Final_con)
    #final_CSV = pd.concat(FINAL_CSV)

    final_con.to_csv(final_csv_path, index=False)

    #features_file for users
    script_dir = sys.path[0]
    results_folder_path = os.path.join(script_dir, "Analysis_Results")

    #Physiochemical parameters
    physico_filename = 'physiochemical_parameters.csv'
    column_names_iv = ['Protein_ID', 'aa_no', 'pH_charge', 'isoelectric point', 'Hydropathy_gravy', 'instability_index', 'aromaticity', '_SecondaryStrC2', '_PolarityC2', '_PolarizabilityT12', '_SecondaryStrT13', '_PolarizabilityD1001', '_PolarizabilityD1025', '_PolarizabilityD1100', '_PolarizabilityD2001', '_PolarizabilityD2100', '_PolarizabilityD3001', '_PolarizabilityD3025', '_PolarizabilityD3100', '_SolventAccessibilityD1025', '_SolventAccessibilityD1100', '_SolventAccessibilityD2001', '_SolventAccessibilityD2025', '_SolventAccessibilityD2100', '_SolventAccessibilityD3001', '_SolventAccessibilityD3025', '_SolventAccessibilityD3100', '_SecondaryStrD1100', '_SecondaryStrD2025', '_SecondaryStrD2100', '_SecondaryStrD3001', '_SecondaryStrD3100', '_ChargeD1001', '_ChargeD1025', '_ChargeD2100', '_ChargeD3100', '_PolarityD1025', '_PolarityD2100', '_PolarityD3100', 'GearyAuto_Hydrophobicity3', 'GearyAuto_Hydrophobicity4', 'GearyAuto_Hydrophobicity7', 'GearyAuto_Hydrophobicity9', 'GearyAuto_Hydrophobicity10', 'GearyAuto_Hydrophobicity14', 'GearyAuto_Hydrophobicity17', 'GearyAuto_Hydrophobicity24', 'GearyAuto_Hydrophobicity25', 'GearyAuto_Hydrophobicity27', 'GearyAuto_Hydrophobicity29', 'GearyAuto_Hydrophobicity30', 'GearyAuto_Polarizability1', 'GearyAuto_Polarizability3', 'GearyAuto_Polarizability5', 'GearyAuto_Polarizability7', 'GearyAuto_Polarizability9', 'GearyAuto_Polarizability10', 'GearyAuto_Polarizability11', 'GearyAuto_Polarizability12', 'GearyAuto_Polarizability14', 'GearyAuto_Polarizability15', 'GearyAuto_Polarizability17', 'GearyAuto_Polarizability19', 'GearyAuto_Polarizability20', 'GearyAuto_Polarizability21', 'GearyAuto_Polarizability22', 'GearyAuto_Polarizability23', 'GearyAuto_Polarizability25', 'GearyAuto_Polarizability26', 'GearyAuto_Polarizability27', 'GearyAuto_Polarizability29', 'GearyAuto_Steric1', 'GearyAuto_Steric2', 'GearyAuto_Steric3', 'GearyAuto_Steric5', 'GearyAuto_Steric7', 'GearyAuto_Steric9', 'GearyAuto_Steric11', 'GearyAuto_Steric13', 'GearyAuto_Steric15', 'GearyAuto_Steric16', 'GearyAuto_Steric17', 'GearyAuto_Steric19', 'GearyAuto_Steric21', 'GearyAuto_Steric23', 'GearyAuto_Steric25', 'GearyAuto_Steric26', 'GearyAuto_Steric27', 'GearyAuto_Steric29', 'GearyAuto_Steric30', 'GearyAuto_Mutability1', 'GearyAuto_Mutability2', 'GearyAuto_Mutability4', 'GearyAuto_Mutability8', 'GearyAuto_Mutability11', 'GearyAuto_Mutability16', 'GearyAuto_Mutability20', 'GearyAuto_Mutability23', 'GearyAuto_Mutability25', 'GearyAuto_Mutability27', 'GearyAuto_Mutability28', 'GearyAuto_Mutability29']
    phsico_params_df = pd.read_csv(final_csv_path)
    phsico_params=pd.DataFrame(phsico_params_df,index=phsico_params_df.index,columns=phsico_params_df.columns)
    phsico_params_columns = phsico_params[column_names_iv]
    physico_result_path = os.path.join(results_folder_path, physico_filename)
    phsico_params_columns.to_csv(physico_result_path, index=False)

    #Signal peptides
    signalp_filename = 'signalp_analysis.csv'
    column_names_v = ['Protein_ID', 'signal_peptide_SP', 'signal_peptide_LIPO', 'signal_peptide_TATLIPO', 'signal_peptide_TAT', 'signal_peptide_PILIN', 'signal_peptide_OTHERS']
    signalp_data = [ID, SP, LIPO, TATLIPO, TAT, PILIN, OTHER]
    df_signalp = pd.DataFrame(data=[signalp_data], columns=column_names_v)
    signalp_result_path = os.path.join(results_folder_path, signalp_filename)
    df_signalp.to_csv(signalp_result_path, index=False)


    #RESACLING THE DATA IN FINAL CSV
    rescale_df = pd.read_csv(final_csv_path)
    Protein_ID = rescale_df["Protein_ID"]
    rescale_df = rescale_df.drop(columns=["level_0", "Protein_ID", "index.1"])

    scaler=per. MinMaxScaler(feature_range=(-1, 1))
    rescaledData=scaler.fit_transform(rescale_df)
    rescaledData=pd.DataFrame(rescaledData,index=rescale_df.index,columns=rescale_df.columns)
    rescaledData.insert(0, 'Protein_ID', Protein_ID)
    columns_to_keep = ['Protein_ID', 'aa_no', 'pH_charge', 'isoelectric point', 'Hydropathy_gravy', 'instability_index', 'aromaticity', 'antigenicity_1', 'b_cells_probability_score', 'mhci_probability_score', 'mhci_rank', 'mhcii_rank', 'mhcii_score', 'surface_probability', 'signal_peptide_SP', 'signal_peptide_LIPO', 'signal_peptide_TAT', 'signal_peptide_TATLIPO', 'signal_peptide_PILIN', 'signal_peptide_OTHER', '_SecondaryStrC2', '_PolarityC2', '_PolarizabilityT12', '_SecondaryStrT13', '_PolarizabilityD1001', '_PolarizabilityD1025', '_PolarizabilityD1100', '_PolarizabilityD2001', '_PolarizabilityD2100', '_PolarizabilityD3001', '_PolarizabilityD3025', '_PolarizabilityD3100', '_SolventAccessibilityD1025', '_SolventAccessibilityD1100', '_SolventAccessibilityD2001', '_SolventAccessibilityD2025', '_SolventAccessibilityD2100', '_SolventAccessibilityD3001', '_SolventAccessibilityD3025', '_SolventAccessibilityD3100', '_SecondaryStrD1100', '_SecondaryStrD2025', '_SecondaryStrD2100', '_SecondaryStrD3001', '_SecondaryStrD3100', '_ChargeD1001', '_ChargeD1025', '_ChargeD2100', '_ChargeD3100', '_PolarityD1025', '_PolarityD2100', '_PolarityD3100', 'GearyAuto_Hydrophobicity3', 'GearyAuto_Hydrophobicity4', 'GearyAuto_Hydrophobicity7', 'GearyAuto_Hydrophobicity9', 'GearyAuto_Hydrophobicity10', 'GearyAuto_Hydrophobicity14', 'GearyAuto_Hydrophobicity17', 'GearyAuto_Hydrophobicity24', 'GearyAuto_Hydrophobicity25', 'GearyAuto_Hydrophobicity27', 'GearyAuto_Hydrophobicity29', 'GearyAuto_Hydrophobicity30', 'GearyAuto_Polarizability1', 'GearyAuto_Polarizability3', 'GearyAuto_Polarizability5', 'GearyAuto_Polarizability7', 'GearyAuto_Polarizability9', 'GearyAuto_Polarizability10', 'GearyAuto_Polarizability11', 'GearyAuto_Polarizability12', 'GearyAuto_Polarizability14', 'GearyAuto_Polarizability15', 'GearyAuto_Polarizability17', 'GearyAuto_Polarizability19', 'GearyAuto_Polarizability20', 'GearyAuto_Polarizability21', 'GearyAuto_Polarizability22', 'GearyAuto_Polarizability23', 'GearyAuto_Polarizability25', 'GearyAuto_Polarizability26', 'GearyAuto_Polarizability27', 'GearyAuto_Polarizability29', 'GearyAuto_Steric1', 'GearyAuto_Steric2', 'GearyAuto_Steric3', 'GearyAuto_Steric5', 'GearyAuto_Steric7', 'GearyAuto_Steric9', 'GearyAuto_Steric11', 'GearyAuto_Steric13', 'GearyAuto_Steric15', 'GearyAuto_Steric16', 'GearyAuto_Steric17', 'GearyAuto_Steric19', 'GearyAuto_Steric21', 'GearyAuto_Steric23', 'GearyAuto_Steric25', 'GearyAuto_Steric26', 'GearyAuto_Steric27', 'GearyAuto_Steric29', 'GearyAuto_Steric30', 'GearyAuto_Mutability1', 'GearyAuto_Mutability2', 'GearyAuto_Mutability4', 'GearyAuto_Mutability8', 'GearyAuto_Mutability11', 'GearyAuto_Mutability16', 'GearyAuto_Mutability20', 'GearyAuto_Mutability23', 'GearyAuto_Mutability25', 'GearyAuto_Mutability27', 'GearyAuto_Mutability28', 'GearyAuto_Mutability29']
    final_df = rescaledData[columns_to_keep]

    # concatenate the script directory with the filename to create the full path of the CSV file
    scaled_csv_path = os.path.join(script_dir, "Rescaled_CSV.csv")
    final_df.to_csv(scaled_csv_path, index=False)
    completed_steps += 1
    progress=int((completed_steps / total_steps) * 100)
    progress_callback("Step 8: Dataset Created Successfully", 8)


def upload_sequence(request):
    if request.method == "POST":
        sequence = request.POST.get("sequence")
        #sequence = sequence.replace("\n", "")
        file = request.FILES.get("file")

        try:
            if sequence:
                script_dir = os.path.dirname(os.path.abspath(__file__))
                file_path = os.path.join(script_dir, "sequences.fasta")

                with open(file_path, "w") as f:
                    f.write(sequence)
            elif file:
                script_dir = os.path.dirname(os.path.abspath(__file__))
                file_path = os.path.join(script_dir, "sequences.fasta")
                with open(file_path, "wb") as f:
                    for chunk in file.chunks():
                        f.write(chunk)

            calculate_features(file_path, progress_callback)
            

        except Exception as e:
            print(f"Error analyzing sequence: {e}")
            return JsonResponse({"status": "Error analyzing sequence"}, status=500)

        return redirect("results")
    
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    log_file_path = os.path.join(script_dir, "vacsol_ai.log")
    file_handler = logging.FileHandler(log_file_path)
    file_handler.setLevel(logging.INFO)
    formatter = ProgressBarFormatter('%(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    with open(log_file_path, 'r') as log_file:
        logs = log_file.readlines()

    # Clear the log file after reading the messages
    with open(log_file_path, 'w') as log_file:
        log_file.truncate(0)

    return render(request, 'upload.html', {'logs': logs})

class ProgressBarFormatter(logging.Formatter):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.completed_steps = 0
        self.total_steps = 10
        self.progress_bar_added = True

    def format(self, record):
        # Check if the log message contains progress information
        if hasattr(record, 'progress'):
            self.completed_steps = record.progress

        # Calculate progress percentage
        progress = int((self.completed_steps / self.total_steps) * 100)

        # Create the progress bar string
        progress_bar = f"[{'#' * (progress // 10)}{'-' * ((100 - progress) // 10)}] {progress}%"

        # Check if the progress bar has already been added
        if not self.progress_bar_added:
            # Modify the log message to include the progress bar
            record.msg = f"{progress_bar}\n{record.msg}"
            self.progress_bar_added = True
        else:
            # Remove the progress bar from subsequent log messages
            record.msg = f"\n{record.msg}"

        return super().format(record)

def get_results(request):
    try:
        script_dir = sys.path[0]
        model_path = os.path.join(script_dir, "model.pkl")

        with open(model_path, "rb") as file:
            model = pickle.load(file)

        csv_path = os.path.join(script_dir, "Rescaled_CSV.csv")
        df = pd.read_csv(csv_path)
        feature_cols = df.drop(labels=["Protein_ID"], axis=1)

        x_new = feature_cols

        out = model.predict(x_new)
        results = pd.DataFrame({"Protein_ID": df["Protein_ID"], "Prediction": out})

        results_data = [
            {"Protein_ID": row["Protein_ID"], "Prediction": int(row["Prediction"])}
            for _, row in results.iterrows()
        ]

        return render(request, "results.html", {"results": results_data})

    except Exception as e:
        print(f"Error analyzing sequence: {e}")
        return JsonResponse({"status": "Error analyzing sequence"}, status=500)
    
def download_csv(request):
    if request.method == 'POST':
        selected_files = request.POST.getlist('csv_files')
        if selected_files:
            # Set the response content type as a zip file
            response = HttpResponse(content_type='application/zip')
            response['Content-Disposition'] = 'attachment; filename="analysis_results.zip"'
            # Create a zip file to store the selected CSV files
            with zipfile.ZipFile(response, 'w') as zip_file:
                # Path to the directory containing the CSV files
                script_dir = sys.path[0]
                directory = os.path.join(script_dir, 'analysis_results')
                for file_name in selected_files:
                    # Path to the CSV file
                    file_path = os.path.join(directory, file_name)
                    # Add the CSV file to the zip file
                    zip_file.write(file_path, os.path.basename(file_path))
            return response
    else:
        return render(request, 'results.html')
    
def delete_files(request):
    script_dir = sys.path[0]
    folder_path = os.path.join(script_dir, 'biolib_results')

    # Delete all files in the folder
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if os.path.isfile(file_path):
            os.remove(file_path)

    return HttpResponseRedirect('/')


### VacSol_ai ###

def sequence_upload(request):
    if request.method == "POST":
        sequence = request.POST.get("sequence")
        #sequence = sequence.replace("\n", "")
        file = request.FILES.get("file")

        try:
            if sequence:
                script_dir = os.path.dirname(os.path.abspath(__file__))
                file_path = os.path.join(script_dir, "sequences.fasta")

                with open(file_path, "w") as f:
                    f.write(sequence)
            elif file:
                script_dir = os.path.dirname(os.path.abspath(__file__))
                file_path = os.path.join(script_dir, "sequences.fasta")
                with open(file_path, "wb") as f:
                    for chunk in file.chunks():
                        f.write(chunk)

            features_annotate(file_path, progress_callback)

        except Exception as e:
            print(f"Error analyzing sequence: {e}")
            return JsonResponse({"status": "Error analyzing sequence"}, status=500)

        return redirect("predict_results")
    
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    log_file_path = os.path.join(script_dir, "vacsol_ai.log")
    file_handler = logging.FileHandler(log_file_path)
    file_handler.setLevel(logging.INFO)
    formatter = ProgressBarFormatter('%(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    with open(log_file_path, 'r') as log_file:
        logs = log_file.readlines()

    # Clear the log file after reading the messages
    with open(log_file_path, 'w') as log_file:
        log_file.truncate(0)

    return render(request, 'upload_sequence.html', {'logs': logs})

def features_annotate(file_path, progress_callback):
    total_steps = 10
    completed_steps = 0
    # Use os.path.abspath to get the absolute file path
    script_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(script_dir, "sequences.fasta")
    file_path = os.path.abspath(file_path)
    MERGED_DATAFRAMES = []
    with open(file_path) as fasta_file:
        parser = fastaparser.Reader(fasta_file)
        seq_num = 1
        ID_list = []
        list_data = []
        for seq in parser:
            ID = seq.id
            sequence = seq.formatted_fasta()
            sequence1 = seq.sequence_as_string()
            ID_list.append(ID)
            df_ID = pd.DataFrame(ID_list, columns=['ID'])
            completed_steps += 1
            progress=int((completed_steps / total_steps) * 100)
            progress_callback(f"Step 1: Sequence {ID} uploaded successfully", 1)

            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type AAC --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature1.tsv')
            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type GAAC --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature2.tsv')
            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type Moran --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature3.tsv')
            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type Geary --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature4.tsv')
            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type NMBroto --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature5.tsv')
            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type CTDC --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature6.tsv')
            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type CTDD --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature7.tsv')
            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type CTDT --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature8.tsv')
            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type CTriad --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature9.tsv')
            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type KSCTriad --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature10.tsv')
            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type SOCNumber --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature11.tsv')
            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type QSOrder --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature12.tsv')
            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type PAAC --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature13.tsv')
            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type APAAC --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature14.tsv')
            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type CKSAAP --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature15.tsv')
            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type DPC --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature16.tsv')
            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type DDE --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature17.tsv')
            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type TPC --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature18.tsv')
            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type CKSAAGP --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature19.tsv')
            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type GDPC --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature20.tsv')
            subprocess.run('python E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/iFeature/iFeature.py --file E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/sequences.fasta --type GTPC --out E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature21.tsv')
            
            completed_steps += 1
            progress=int((completed_steps / total_steps) * 100)
            progress_callback("Step 2: Physicochemical analysis completed", 2)


            df1 = pd.read_table('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature1.tsv')
            df2 = pd.read_table('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature2.tsv')
            df3 = pd.read_table ('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature3.tsv')
            df4 = pd.read_table('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature4.tsv')
            df5 = pd.read_table('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature5.tsv')
            df6 = pd.read_table ('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature6.tsv')
            df7 = pd.read_table('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature7.tsv')
            df8 = pd.read_table('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature8.tsv')
            df9 = pd.read_table ('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature9.tsv')
            df10 = pd.read_table('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature10.tsv')
            df11 = pd.read_table('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature11.tsv')
            df12 = pd.read_table ('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature12.tsv')
            df13 = pd.read_table('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature13.tsv')
            df14 = pd.read_table ('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature14.tsv')
            df15 = pd.read_table('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature15.tsv')
            df16 = pd.read_table ('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature16.tsv')
            df17 = pd.read_table('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature17.tsv')
            df18 = pd.read_table('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature18.tsv')
            df19 = pd.read_table ('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature19.tsv')
            df20 = pd.read_table('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature20.tsv')
            df21 = pd.read_table ('E:/ASAB/VacSol-AI/ESKAPE/Real_time_server/upload/ifeature21.tsv')

            MERGED_DATAFRAMES = []
            df_all = [df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15, df16, df17, df18, df19, df20, df21]
            con1 = pd.concat(df_all, axis="columns")
            con1 = con1.loc[:, ~con1.columns.duplicated()] 
            MERGED_DATAFRAMES.append(con1)

            completed_steps += 1
            progress=int((completed_steps / total_steps) * 100)
            progress_callback("Step 3: Creating Dataset", 3)

            con1.to_csv ('testifeature.csv', index=False)

            df_final = pd.read_csv('testifeature.csv')
            columns_to_keep = ['A', 'CE.gap0', 'EY.gap0', 'GD.gap0', 'IE.gap0', 'LD.gap0', 'LI.gap0', 'NM.gap0', 'PK.gap0', 'TE.gap0', 'YP.gap0', 'MG.gap1', 'NM.gap1', 'RG.gap1', 'TV.gap1', 'TW.gap1', 'YK.gap1', 'AF.gap2', 'EI.gap2', 'FN.gap2', 'FR.gap2', 'GF.gap2', 'MD.gap2', 'NA.gap2', 'SM.gap2', 'SW.gap2', 'TQ.gap2', 'YE.gap2', 'AI.gap3', 'AY.gap3', 'IG.gap3', 'IH.gap3', 'PF.gap3', 'RF.gap3', 'VS.gap3', 'FV.gap4', 'KY.gap4', 'LS.gap4', 'NA.gap4', 'PL.gap4', 'PS.gap4', 'YH.gap4', 'DM.gap5', 'FP.gap5', 'GH.gap5', 'GR.gap5', 'HI.gap5', 'HM.gap5', 'HN.gap5', 'PN.gap5', 'SF.gap5', 'AAH', 'AET', 'AHQ', 'AIT', 'AMV', 'AQE', 'AVQ', 'AVT', 'AYA', 'DDN', 'DKH', 'DNR', 'EID', 'ELA', 'EMS', 'EQT', 'ESK', 'ETK', 'EVG', 'FDH', 'FIK', 'GDE', 'GGI', 'GGL', 'GLL', 'GWE', 'HVK', 'IAI', 'IIT', 'ILT', 'IQF', 'IVN', 'IYS', 'KGS', 'KIN', 'KKL', 'KVN', 'LDF', 'LDS', 'LKN', 'LNG', 'LNY', 'LSD', 'NAA', 'NGK', 'NNF', 'NYK', 'PLA', 'PVK', 'QND', 'RID', 'SAI', 'SKD', 'SRM', 'SSA', 'STL', 'TEA', 'TKD', 'TMM', 'TNE', 'TSG', 'TTV', 'TYK', 'VEK', 'VNG', 'VRN', 'VSA', 'VYN', 'YAN', 'YGA', 'YGY', 'YSP', 'CIDH920105.lag1', 'CIDH920105.lag26', 'CHAM820102.lag3', 'CHAM820102.lag15', 'CHAM810101.lag14', 'CHAM810101.lag24', 'g1.g1.g3', 'g1.g4.g4', 'g1.g4.g5', 'g1.g6.g4', 'g3.g5.g2', 'g4.g5.g3']
            final_df = df_final[columns_to_keep]
            final_df.to_csv('ifeatures.csv', index=False)

            #features_file for users
            script_dir = sys.path[0]
            results_folder_path = os.path.join(script_dir, "Analysis_Results")

            #physicochemical features
            physico_filename = 'physicochemical_parameters.csv'
            physico_result_path = os.path.join(results_folder_path, physico_filename)
            final_df.to_csv(physico_result_path, index=False)


            sprob_res = iedb.query_bcell_epitope(method="Emini", sequence= sequence1, window_size=9)
            sprob_res['Score'] = pd.to_numeric(sprob_res['Score'], errors='coerce')
            sprob = sprob_res["Score"].mean()

            completed_steps += 1
            progress=int((completed_steps / total_steps) * 100)
            progress_callback("Step 4: Adhesion Probability analysis completed", 4)

            add_scores2 = pd.read_csv('ifeatures.csv')
            add_scores2['SPAAN-Pad-value'] = sprob
            add_scores2.insert(0, 'Protein_ID', df_ID)

            completed_steps += 1
            progress=int((completed_steps / total_steps) * 100)
            progress_callback("Step 5: Added Adhesion Probability scores to Dataset", 5)

            #adhesion_probability feature
            adhesion_filename = 'adhesion_probability.csv'
            column_names_iii = ['Protein_ID' , 'adhesion_probability']
            adhesion_data = [ID, sprob]
            df = pd.DataFrame(data=[adhesion_data], columns=column_names_iii)
            df_adhesion = pd.DataFrame(data=[adhesion_data], columns=column_names_iii)
            adhesion_result_path = os.path.join(results_folder_path, adhesion_filename)
            df_adhesion.to_csv(adhesion_result_path, index=False)


            add_scores2.to_csv('ifeatures_updated.csv', index=False)

    signalp_6 = biolib.load('DTU/SignalP_6')
    signalp6_job = signalp_6.cli(args='--fastafile vacsol_ai/sequences.fasta --output_dir output')
    signalp6_job.save_files('biolib_results/signalP')

    script_dir = sys.path[0]
    signalP_folder_path = os.path.join(script_dir, "biolib_results/signalP")

    completed_steps += 1
    progress=int((completed_steps / total_steps) * 100)
    progress_callback("Step 6: Signal Peptide Analysis Completed", 6)

    #read signalp results:    
    sp_table_path = os.path.join(signalP_folder_path, "prediction_results.txt")
    df = pd.read_table(sp_table_path, sep="\t", header=None, skiprows=1)
    df.columns = ["ID", "Prediction", "OTHER", "SP(Sec/SPI)", "LIPO(Sec/SPII)", "TAT(Tat/SPI)", "TATLIPO(Sec/SPII)", "PILIN(Sec/SPIII)", "CS Position"]
    df.drop(index=df.index[0], axis=1, inplace=True)

    #extract the scores
    SP = df['SP(Sec/SPI)'].values
    LIPO = df['LIPO(Sec/SPII)'].values
    TAT = df['TAT(Tat/SPI)'].values
    TATLIPO = df['TATLIPO(Sec/SPII)'].values
    PILIN = df['PILIN(Sec/SPIII)'].values
    OTHER = df['OTHER'].values

    #features_file for users
    script_dir = sys.path[0]
    results_folder_path = os.path.join(script_dir, "Analysis_Results")

    #Signal peptides features
    signalp_filename = 'signalp_analysis.csv'
    column_names_v = ['Protein_ID', 'signal_peptide_SP', 'signal_peptide_LIPO', 'signal_peptide_TATLIPO', 'signal_peptide_TAT', 'signal_peptide_PILIN', 'signal_peptide_OTHERS']
    signalp_data = [ID, SP, LIPO, TATLIPO, TAT, PILIN, OTHER]
    df_signalp = pd.DataFrame(data=[signalp_data], columns=column_names_v)
    signalp_result_path = os.path.join(results_folder_path, signalp_filename)
    df_signalp.to_csv(signalp_result_path, index=False)


    deeptmhmm = biolib.load('DTU/DeepTMHMM')
    deeptmhmm_job = deeptmhmm.cli(args='--fasta vacsol_ai/sequences.fasta')
    deeptmhmm_job.save_files('biolib_results/deeptmhmm')


    completed_steps += 1
    progress=int((completed_steps / total_steps) * 100)
    progress_callback("Step 7: Transmembrane Helix prediction completed", 7)

    directory = 'biolib_results/deeptmhmm'

    files = os.listdir(directory)

    csv_files = [file for file in files if file.endswith('.csv')]

    csv_file_path = os.path.join(directory, csv_files[0])
    df = pd.read_csv(csv_file_path, skiprows=[0])
    df = df.iloc[1:].reset_index(drop=True)

    df['Beta'] = pd.to_numeric(df['Beta'], errors='coerce')
    TMRs_TMB = df['Beta'].mean()

    df['Membrane'] = pd.to_numeric(df['Membrane'], errors='coerce')
    TMRs_TMH = df['Membrane'].mean()

    df['Signal'] = pd.to_numeric(df['Signal'], errors='coerce')
    SignalP = df['Signal'].mean()

    #DeepTMHMM
    deeptmhmm_filename = 'deeptmhmm.csv'
    column_names_v = ['Protein_ID', 'DeepTMHMM-No. of TMRs-TMB', 'DeepTMHMM-No. of TMRs-TMH', 'DeepTMHMM-SP']
    deeptmhmm_data = [ID, TMRs_TMB, TMRs_TMH, SignalP]
    df_deeptmhmm = pd.DataFrame(data=[deeptmhmm_data], columns=column_names_v)
    deeptmhmm_result_path = os.path.join(results_folder_path, deeptmhmm_filename)
    df_deeptmhmm.to_csv(deeptmhmm_result_path, index=False)


    completed_steps += 1
    progress=int((completed_steps / total_steps) * 100)
    progress_callback("Step 8: Adding SignalP Scores to Dataset", 8)

    add_scores = pd.read_csv('ifeatures_updated.csv ')
    add_scores["signal_peptide_SP"] = SP
    add_scores["signal_peptide_LIPO"] = LIPO
    add_scores["signal_peptide_TAT"] = TAT
    add_scores["signal_peptide_TATLIPO"] = TATLIPO
    add_scores["signal_peptide_PILIN"] = PILIN
    add_scores["signal_peptide_OTHER"] = OTHER

    completed_steps += 1
    progress=int((completed_steps / total_steps) * 100)
    progress_callback("Step 9: Adding Transmembrane Helix Scores to Dataset", 9)

    add_scores["DeepTMHMM-No. of TMRs-TMB"] = TMRs_TMB
    add_scores["DeepTMHMM-No. of TMRs-TMH"] = TMRs_TMH
    add_scores["DeepTMHMM-SP"] = SignalP

    add_scores.to_csv('finalifeatures_updated.csv', index=False)

    completed_steps += 1
    progress=int((completed_steps / total_steps) * 100)
    progress_callback("Step 10: Dataset Creating Successfully, Starting AI Predictions", 10)


def predict_results(request):
    try:
        script_dir = sys.path[0]
        model_path = os.path.join(script_dir, "BACmodel.pkl")

        with open(model_path, "rb") as file:
            model = pickle.load(file)

        csv_path = os.path.join(script_dir, "finalifeatures_updated.csv")
        df = pd.read_csv(csv_path)
        feature_cols = df.drop(labels=["Protein_ID"], axis=1)

        x_new = feature_cols

        out = model.predict(x_new)
        results = pd.DataFrame({"Protein_ID": df["Protein_ID"], "Prediction": out})

        results_data = [
            {"Protein_ID": row["Protein_ID"], "Prediction": int(row["Prediction"])}
            for _, row in results.iterrows()
        ]

        return render(request, "results_vacsolai.html", {"results": results_data})

    except Exception as e:
        print(f"Error analyzing sequence: {e}")
        return JsonResponse({"status": "Error analyzing sequence"}, status=500)
    
def delete_files_vacsolai(request):
    script_dir = sys.path[0]
    folder_path = os.path.join(script_dir, 'biolib_results')

    # Delete all files in the folder
    for root, dirs, files in os.walk(folder_path, topdown=False):
        for file_name in files:
            file_path = os.path.join(root, file_name)
            os.remove(file_path)

        for dir_name in dirs:
            dir_path = os.path.join(root, dir_name)
            shutil.rmtree(dir_path)

    return HttpResponseRedirect('/')