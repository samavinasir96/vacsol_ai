import pickle
from django.conf import settings
import iedb
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from PyBioMed.PyPretreat import PyPretreatPro
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from PyBioMed import Pyprotein
import pandas as pd
import os
import fastaparser
from pandas import DataFrame
from sklearn import preprocessing as per
import sys    
from django.http import HttpResponse, JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.views import View
from django.shortcuts import redirect, render
import pickle
import pandas as pd
import os
import sys


@csrf_exempt
def home(request):
    return render (request, 'index.html')

def progress_callback(message):
    sys.stdout.write(message + '\n')
    sys.stdout.flush()

def calculate_features(file_path, progress_callback, run_analysis=False):
    # Use os.path.abspath to get the absolute file path
    file_path = os.path.abspath(file_path)
    LIST_RESULTS_1 = []
    LIST_RESULTS_2 = []
    LIST_RESULTS_3 = []
    LIST_RESULTS_4 = []
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

            # Send POST request to MHC class I peptide binding prediction tool:
            progress_callback(f"Sending sequence {ID} to Epitope Analysis Predictions...")
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

            # Send POST request to surface probability prediction tool:
            progress_callback(f"Sending sequence {ID} to Adhesion Probability Predictions...")
            sprob_res = iedb.query_bcell_epitope(method="Emini", sequence= sequence1, window_size=9)

            # Send POST request to antigenicity prediction tool:
            progress_callback(f"Sending sequence {ID} to Antigenicity Predictions...")
            antigenicity_1 = iedb.query_bcell_epitope(method="Kolaskar-Tongaonkar", sequence= sequence1, window_size=9)

            # Getting means - mhci
            mhci_res['score'] = pd.to_numeric(mhci_res['score'], errors='coerce')
            mhci_res['percentile_rank'] = pd.to_numeric(mhci_res['percentile_rank'], errors='coerce')

            df_score_mhci = mhci_res["score"].mean()
            df_rank_mhci = mhci_res["percentile_rank"].mean()

            # Getting means - mhcii
            mhcii_res['score'] = pd.to_numeric(mhcii_res['score'], errors='coerce')
            mhcii_res['rank'] = pd.to_numeric(mhcii_res['rank'], errors='coerce')
#            mhcii_res['adjusted_rank'] = pd.to_numeric(mhcii_res['adjusted_rank'], errors='coerce')
 #           mhcii_res['comblib_score'] = pd.to_numeric(mhcii_res['comblib_score'], errors='coerce')
  #          mhcii_res['comblib_adjusted_rank'] = pd.to_numeric(mhcii_res['comblib_adjusted_rank'], errors='coerce')
   #         mhcii_res['smm_align_adjusted_rank'] = pd.to_numeric(mhcii_res['smm_align_adjusted_rank'], errors='coerce')
  #          mhcii_res['nn_align_ic50'] = pd.to_numeric(mhcii_res['nn_align_ic50'], errors='coerce')
   #         mhcii_res['nn_align_adjusted_rank'] = pd.to_numeric(mhcii_res['nn_align_adjusted_rank'], errors='coerce')

            rank_mhcii = mhcii_res["rank"].mean()
            score_mhcii = mhcii_res["score"].mean()
 #           comblib_score_mhcii = mhcii_res['comblib_score'].mean()
 #           comblib_rank_mhcii = mhcii_res["comblib_adjusted_rank"].mean()
  #          smm_align_ic50_mhcii = mhcii_res["smm_align_ic50"].mean()
   #         smm_align_rank_mhcii = mhcii_res["smm_align_adjusted_rank"].mean()
    #        nn_align_ic50_mhcii = mhcii_res["nn_align_ic50"].mean()
     #       nn_align_rank_mhcii = mhcii_res["nn_align_adjusted_rank"].mean()

            # Getting means - bcells
            bcell_res['Score'] = pd.to_numeric(bcell_res['Score'], errors='coerce')
            df_bcells_final = bcell_res["Score"].mean()

            # Getting means - s_probability
            sprob_res['Score'] = pd.to_numeric(sprob_res['Score'], errors='coerce')
            df_sprob_final = sprob_res["Score"].mean()

            # Getting means - antigenicity(scale1)
            antigenicity_1['Score'] = pd.to_numeric(antigenicity_1['Score'], errors='coerce')
            antigenicity_1_final = antigenicity_1["Score"].mean()
            
            #Analysis of Physiochemical Features:
            progress_callback(f"Sending sequence {ID} to Physicochemical Parameters Predictions...")
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
            progress_callback(f"Sending sequence {ID} to caculate feature descriptors...")
            aa_no= PyPretreatPro.ProteinCheck(sequence1)
            protein_class = Pyprotein.PyProtein(sequence1)
            descriptor1 = protein_class.GetCTD()
            descriptor2 = protein_class.GetGearyAuto()
            descriptor3 = protein_class.GetMoranAuto()
            descriptor4 = protein_class.GetMoreauBrotoAuto()
            
            #df1 = pd.DataFrame.from_dict(descriptor1)
            df1 = pd.DataFrame.from_dict(descriptor1,orient='index').transpose()
            df2 = pd.DataFrame.from_dict(descriptor2,orient='index').transpose()
            df3 = pd.DataFrame.from_dict(descriptor3,orient='index').transpose()
            df4 = pd.DataFrame.from_dict(descriptor4,orient='index').transpose()
            
            df_all = [df1, df2, df3, df4]
            con1 = pd.concat(df_all, axis="columns")
            MERGED_DATAFRAMES.append(con1)

            #Raw Results files
            progress_callback(f"Creating DataSet...")
            """
            #i: Physiochemical csv
            column_names_i = ['Protein_ID','aa_no','molecular_weight', 'pH_charge', 'isoelectric point', 'Hydropathy_gravy', 'instability_index', 'aromaticity']
            physiochemical_data = [ID, aa_no, scale1, scale2, scale3, scale4, scale5, scale6]
            df_physiochemical = pd.DataFrame(data=[physiochemical_data], columns=column_names_i)
            df_physiochemical.to_csv('E:/ASAB/VacSol-AI/VacSol-ML-ESKAPE-/physiochemical_analysis.csv', index=False)

            #ii: Epitope_analysis csv
            column_names_ii = ['Protein_ID', 'b_cells_probability_score', 'mhci_probability_score', 'mhci_rank', 'mhcii_rank', 'mhcii_comblib_score', 'mhcii_comblib_rank', 'mhcii_smm_align_ic50', 'mhcii_smm_align_rank', 'mhcii_nn_align_ic50', 'mhcii_nn_align_rank']
            epitope_data = [ID, df_bcells_final, df_score_mhci, df_rank_mhci, rank_mhcii, comblib_score_mhcii, comblib_rank_mhcii, smm_align_ic50_mhcii, smm_align_rank_mhcii, nn_align_ic50_mhcii, nn_align_rank_mhcii]
            df_epitope = pd.DataFrame(data=[epitope_data], columns=column_names_ii)
            df_epitope.to_csv('E:/ASAB/VacSol-AI/VacSol-ML-ESKAPE-/epitope_analysis.csv', index=False)

            #iii: adhesion_probability
            column_names_iii = ['Protein_ID', 'surface_probability']
            adhesion_data = [ID, df_sprob_final]
            df = pd.DataFrame(data=[adhesion_data], columns=column_names_iii)
            df_adhesion = pd.DataFrame(data=[adhesion_data], columns=column_names_iii)
            df_adhesion.to_csv('E:/ASAB/VacSol-AI/VacSol-ML-ESKAPE-/adhesion_analysis.csv', index=False)

            #iv: antigenicity
            column_names_iv = ['Protein_ID','antigenicity_1']
            antigenicity_data = [ID, antigenicity_1_final]
            df = pd.DataFrame(data=[antigenicity_data], columns=column_names_iv)
            df_antigenicity = pd.DataFrame(data=[antigenicity_data], columns=column_names_iv)
            df_antigenicity.to_csv('E:/ASAB/VacSol-AI/VacSol-ML-ESKAPE-/antigenicity_analysis.csv', index=False)
            """

    con = pd.concat(list_data).reset_index()
    print(con)

    # concatenate the script directory with the filename to create the full path of the CSV file
    script_dir = sys.path[0]
    final_csv_path = os.path.join(script_dir, "FINAL.csv")


    con.to_csv(final_csv_path, index=False)

    #prediction of signal peptides (SignalP. version 6.0)
    # get the path of the directory where the script is located
    script_dir = sys.path[0]

    # concatenate the script directory with the filename to create the full path of the FASTA file
    fasta_file_path = os.path.join(script_dir, "sequences.fasta")

    # create the output directory path
    csv_folder_path = os.path.join(script_dir, "Analysis_Results")
    os.makedirs(csv_folder_path, exist_ok=True)

    output_dir_path = os.path.join(csv_folder_path)

    #prediction of signal peptides (SignalP. version 6.0)
    os.system(f'signalp6 --fastafile {fasta_file_path} --organism other --output_dir {output_dir_path}')

    #read signalp results:
    progress_callback(f"Sending sequence {ID} to Signal Peptide Predictions...")
    
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
    
    progress_callback(f"Creating DataSet...")
    add_scores = pd.read_csv(final_csv_path)
    add_scores["signal_peptide_SP"] = SP
    add_scores["signal_peptide_LIPO"] = LIPO
    add_scores["signal_peptide_TAT"] = TAT
    add_scores["signal_peptide_TATLIPO"] = TATLIPO
    add_scores["signal_peptide_PILIN"] = PILIN
    add_scores["signal_peptide_OTHER"] = OTHER

    MERGED = pd.concat(MERGED_DATAFRAMES).reset_index()

    data_CSV = [add_scores, MERGED]
    final_con = pd.concat(data_CSV, axis="columns").reset_index()
    #FINAL_CSV.append(Final_con)
    #final_CSV = pd.concat(FINAL_CSV)

    final_con.to_csv(final_csv_path, index=False)

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
    progress_callback(f"Features calculated successfully!")

def upload_sequence(request):
    if request.method == "POST":
        sequence = request.POST.get("sequence")
        sequence = sequence.replace("\n", "")
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

        except Exception as e:
            print(f"Error analyzing sequence: {e}")
            return JsonResponse({"status": "Error analyzing sequence"}, status=500)

        return redirect("results")

    return render(request, "upload.html")


def get_results(request):
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
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
        
    
def complete_analysis_view(request):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    directory_path = os.path.join(settings.MEDIA_ROOT, script_dir, "Analysis_Results")  # Set the path to your directory

    # Ensure the directory path is valid
    if os.path.exists(directory_path) and os.path.isdir(directory_path):
        # Delete all files in the directory
        for file_name in os.listdir(directory_path):
            file_path = os.path.join(directory_path, file_name)
            if os.path.isfile(file_path):
                os.remove(file_path)

        return HttpResponse("Analysis completed. All files deleted.")
    else:
        return HttpResponse("Invalid directory path.")


