#####################################
#                                   #
#   Script by Marion PERRIER        #
#   marion.perrier@inserm.fr        #
#                                   #
#   Feel free to use, modify        #
#   or redistribute                 #
#                                   #
#####################################

'''
This script is to clean the BCR and TCR data from all_contig_annotation.csv provided by 10XGenomics.
It will create a more compact .csv file :
    - Only full length (full_length == True)
    - Only cells (is_cell == True)
    - One row = one cell.
'''

#!/usr/bin/python3
import pandas
import numpy as np
import re

is_new = True

while(is_new) :
    file_path = input("Path ? ")
    #Open the csv_file
    df = pandas.read_csv(file_path)
    df_out = []

    #Delete full_length == false && is_cell == false
    df = df[df.full_length == True]
    df = df[df.is_cell == True]

    #We take the names of the columns
    cell_id_unique = pandas.unique(df.barcode)
    df_out = pandas.DataFrame(columns = ['clonotype', 'IGH_1_v', 'IGH_1_d', 'IGH_1_j', 'IGH_1_c', 'IGL_1_v', 'IGL_1_d', 'IGL_1_j', 'IGL_1_c', 'IGK_1_v', 'IGK_1_d', 'IGK_1_j', 'IGK_1_c'], index = cell_id_unique)

    #Search in the main Dataframe the first index
    for index, df_row in df_out.iterrows() :
        #Creation of a itermediate Dataframe, containing all rows with the same ID cell
        temp_df = df[df.barcode == index]

        #We take the clonotype (first clonotype on the temporary Dataframe)
        df_out.loc[index, 'clonotype'] = temp_df['raw_clonotype_id'].iloc[0]

        #For every row on this intermediate dataframe, we are filling the empty cells
        compteur_idh = 1
        compteur_idl = 1
        compteur_idk = 1

        #Row by row the intermediate DF
        for line_temp_index, line_temp_df in temp_df.iterrows() :
            #Creation of a list with all the columns name
            nom_colonne = list(df_out.columns)

            #Look at the type (IGH, IGL or IGK)
            if line_temp_df["chain"] == 'IGH':
                #Is there more than 2 IGH?
                if len(temp_df[temp_df['chain'] == 'IGH']) > 1 :

                    #If yes : Then we will fill the first column, before adding as much columns as there is IDH chains (does it make any sense?!)
                    #First, we want to check is the column is not already existing
                    r = re.compile('IGH_'+str(compteur_idh)+'.')

                    #If yes : We need to fill it
                    if any(r.match(nom) for nom in nom_colonne) :
                        df_out.loc[index, 'IGH_'+str(compteur_idh)+'_v'] = line_temp_df['v_gene']
                        df_out.loc[index, 'IGH_'+str(compteur_idh)+'_d'] = line_temp_df['d_gene']
                        df_out.loc[index, 'IGH_'+str(compteur_idh)+'_j'] = line_temp_df['j_gene']
                        df_out.loc[index, 'IGH_'+str(compteur_idh)+'_c'] = line_temp_df['c_gene']
                        compteur_idh += 1

                    #If no : We need to create it
                    else :
                        df_out.insert(len(df_out.columns), 'IGH_'+str(compteur_idh)+'_v', np.nan)
                        df_out.insert(len(df_out.columns), 'IGH_'+str(compteur_idh)+'_d', np.nan)
                        df_out.insert(len(df_out.columns), 'IGH_'+str(compteur_idh)+'_j', np.nan)
                        df_out.insert(len(df_out.columns), 'IGH_'+str(compteur_idh)+'_c', np.nan)

                        #Then add values in it
                        df_out.loc[index, 'IGH_'+str(compteur_idh)+'_v'] = line_temp_df['v_gene']
                        df_out.loc[index, 'IGH_'+str(compteur_idh)+'_d'] = line_temp_df['d_gene']
                        df_out.loc[index, 'IGH_'+str(compteur_idh)+'_j'] = line_temp_df['j_gene']
                        df_out.loc[index, 'IGH_'+str(compteur_idh)+'_c'] = line_temp_df['c_gene']
                        compteur_idh += 1

                #If there is only one IGH, we fill the appropriate column
                else :
                    #We add all the genes
                    df_out.loc[index, 'IGH_1_v'] = line_temp_df['v_gene']
                    df_out.loc[index, 'IGH_1_d'] = line_temp_df['d_gene']
                    df_out.loc[index, 'IGH_1_j'] = line_temp_df['j_gene']
                    df_out.loc[index, 'IGH_1_c'] = line_temp_df['c_gene']
                    compteur_idh += 1


            elif line_temp_df["chain"] == 'IGL' :
                #Is there more than 2 IGL?
                if len(temp_df[temp_df['chain'] == 'IGL']) > 1 :
                    nom_colonne = list(df_out)
                    r = re.compile('IGL_'+str(compteur_idl)+'.')
                    if any(r.match(nom) for nom in nom_colonne) :
                        df_out.loc[index, 'IGL_'+str(compteur_idl)+'_v'] = line_temp_df['v_gene']
                        df_out.loc[index, 'IGL_'+str(compteur_idl)+'_d'] = line_temp_df['d_gene']
                        df_out.loc[index, 'IGL_'+str(compteur_idl)+'_j'] = line_temp_df['j_gene']
                        df_out.loc[index, 'IGL_'+str(compteur_idl)+'_c'] = line_temp_df['c_gene']
                        compteur_idl += 1

                    else : #If no, we need to create the columns
                        df_out.insert(len(df_out.columns), 'IGL_'+str(compteur_idl)+'_v', np.nan)
                        df_out.insert(len(df_out.columns), 'IGL_'+str(compteur_idl)+'_d', np.nan)
                        df_out.insert(len(df_out.columns), 'IGL_'+str(compteur_idl)+'_j', np.nan)
                        df_out.insert(len(df_out.columns), 'IGL_'+str(compteur_idl)+'_c', np.nan)

                        #Then add values in it
                        df_out.loc[index, 'IGL_'+str(compteur_idl)+'_v'] = line_temp_df['v_gene']
                        df_out.loc[index, 'IGL_'+str(compteur_idl)+'_d'] = line_temp_df['d_gene']
                        df_out.loc[index, 'IGL_'+str(compteur_idl)+'_j'] = line_temp_df['j_gene']
                        df_out.loc[index, 'IGL_'+str(compteur_idl)+'_c'] = line_temp_df['c_gene']
                        compteur_idl += 1
                else :
                    #We add all genes
                    df_out.loc[index, 'IGL_'+str(compteur_idl)+'_v'] = line_temp_df['v_gene']
                    df_out.loc[index, 'IGL_'+str(compteur_idl)+'_d'] = line_temp_df['d_gene']
                    df_out.loc[index, 'IGL_'+str(compteur_idl)+'_j'] = line_temp_df['j_gene']
                    df_out.loc[index, 'IGL_'+str(compteur_idl)+'_c'] = line_temp_df['c_gene']
                    compteur_idl += 1

            else :
                #Is there more than 2 IGK?
                if len(temp_df[temp_df['chain'] == 'IGK']) > 1 :
                    nom_colonne = list(df_out)
                    r = re.compile('IGK_'+str(compteur_idk)+'.')
                    if any(r.match(nom) for nom in nom_colonne) :
                        df_out.loc[index, 'IGK_'+str(compteur_idk)+'_v'] = line_temp_df['v_gene']
                        df_out.loc[index, 'IGK_'+str(compteur_idk)+'_d'] = line_temp_df['d_gene']
                        df_out.loc[index, 'IGK_'+str(compteur_idk)+'_j'] = line_temp_df['j_gene']
                        df_out.loc[index, 'IGK_'+str(compteur_idk)+'_c'] = line_temp_df['c_gene']
                        compteur_idk += 1

                    else : #If no, we need to create the columns
                        df_out.insert(len(df_out.columns), 'IGK_'+str(compteur_idk)+'_v', np.nan)
                        df_out.insert(len(df_out.columns), 'IGK_'+str(compteur_idk)+'_d', np.nan)
                        df_out.insert(len(df_out.columns), 'IGK_'+str(compteur_idk)+'_j', np.nan)
                        df_out.insert(len(df_out.columns), 'IGK_'+str(compteur_idk)+'_c', np.nan)

                        #Then add values in it
                        df_out.loc[index, 'IGK_'+str(compteur_idk)+'_v'] = line_temp_df['v_gene']
                        df_out.loc[index, 'IGK_'+str(compteur_idk)+'_d'] = line_temp_df['d_gene']
                        df_out.loc[index, 'IGK_'+str(compteur_idk)+'_j'] = line_temp_df['j_gene']
                        df_out.loc[index, 'IGK_'+str(compteur_idk)+'_c'] = line_temp_df['c_gene']
                        compteur_idk += 1
                else :
                    #We add all genes
                    df_out.loc[index, 'IGK_'+str(compteur_idk)+'_v'] = line_temp_df['v_gene']
                    df_out.loc[index, 'IGK_'+str(compteur_idk)+'_d'] = line_temp_df['d_gene']
                    df_out.loc[index, 'IGK_'+str(compteur_idk)+'_j'] = line_temp_df['j_gene']
                    df_out.loc[index, 'IGK_'+str(compteur_idk)+'_c'] = line_temp_df['c_gene']
                    compteur_idk += 1

        compteur_idh = 1
        compteur_idk = 1
        compteur_idl = 1

    print(df_out)

    ask_for_name = input('Name of the file? ')

    #na_rep : Transform all the "NaN" to "None"
    df_out.to_csv(ask_for_name+'.csv', na_rep='None', sep='\t', encoding='utf-8')

    want_a_new = input('New file? o/N ')
    if want_a_new != 'O' and want_a_new != 'o' :
        is_new = False
