import pandas as pd
import plotly.graph_objects as go
from taxontabletools.taxontable_manipulation import strip_metadata
from taxontabletools.taxontable_manipulation import collect_metadata
from taxontabletools.taxontable_manipulation import add_metadata
from taxontabletools.taxontable_manipulation import reduce_taxontable
import PySimpleGUI as sg
from pathlib import Path

def itt_import(TaXon_table_xlsx, ITT_table_xlsx, path_to_outdirs):

    ## read itt sheet
    itt_df = pd.read_excel(Path(ITT_table_xlsx))
    ## rename columns
    itt_df.columns = itt_df.values.tolist()[0]
    ## select trait data
    itt_df = itt_df[1:587]

    ## collect all present families
    itt_families = itt_df['family'].values.tolist()

    ## create a list of all available trait categories
    trait_categories = itt_df.columns.tolist()[3:]

    ## create a dictionary to store the trait of each family
    family_traits = {}

    for family in itt_df.values.tolist():
        ## collect family name
        family_name = family[1]
        ## store traits
        family_traits[family_name] = family[3:]

    ## load TaxonTable
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx).fillna('')
    TaXon_table_df_metadata = collect_metadata(TaXon_table_df)
    TaXon_table_df = strip_metadata(TaXon_table_df)

    OTU_list = []

    for OTU in TaXon_table_df.values.tolist():
        ## collect family name
        family_name = OTU[4]

        traits_list = []

        larva_habitat_dict = {}

        if family_name in family_traits.keys():

            ## TRAIT 1
            larva_habitat = ['terrestrial', 'aquatic']
            t1 = family_traits[family_name][0:2]
            lst = []
            for i,j in enumerate(t1):
                if t1[i] > 0:
                    lst.append(larva_habitat[i])
            traits_list.append('+'.join(lst))

            ## TRAIT 2
            adult_habitat = ['terrestrial', 'aquatic']
            t1 = family_traits[family_name][2:4]
            lst = []
            for i,j in enumerate(t1):
                if t1[i] > 0:
                    lst.append(adult_habitat[i])
            traits_list.append('+'.join(lst))

            ## TRAIT 3
            larva_diet_types = ['phytophagous', 'zoophagous', 'mycetophagous', 'saprophagous', 'detritophagous', 'coprophagous']
            t1 = family_traits[family_name][4:10]
            lst = []
            for i,j in enumerate(t1):
                if t1[i] > 0:
                    lst.append(larva_diet_types[i])
            traits_list.append('+'.join(lst))

            ## TRAIT 4
            larva_diet_types_phytophageous = ['phyllophagous', 'saproxylic', 'sap-sucking', 'stem', 'flower', 'seed', 'gall-inducing', 'miners', 'roots']
            t1 = family_traits[family_name][10:19]
            lst = []
            for i,j in enumerate(t1):
                if t1[i] > 0:
                    lst.append(larva_diet_types_phytophageous[i])
            traits_list.append('+'.join(lst))

            ## TRAIT 5
            larva_diet_types_zoophagus = ['predator', 'micropredator', 'parasite', 'parasitoid', 'necrophorous']
            t1 = family_traits[family_name][20:25]
            lst = []
            for i,j in enumerate(t1):
                if t1[i] > 0:
                    lst.append(larva_diet_types_zoophagus[i])
            traits_list.append('+'.join(lst))

            OTU_list.append(OTU[0:9] + traits_list + OTU[9:])

        else:
            OTU_list.append(OTU[0:9] + ['']*5 + OTU[9:])

    ## create a new dataframe including the ITT traits
    new_colums = TaXon_table_df.columns.tolist()[0:9] + ['larva_habitat', 'adult_habitat', 'larva_diet_types', 'larva_diet_types_phytophageous', 'larva_diet_types_zoophagus'] + TaXon_table_df.columns.tolist()[9:]

    new_df = pd.DataFrame(OTU_list, columns=new_colums)

    ## add already existing metadata back to the df
    if len(TaXon_table_df_metadata.columns) != 1:
        new_df = add_metadata(new_df, TaXon_table_df_metadata)

    ## save dataframe
    output_pdf = Path(str(path_to_outdirs) + "/TaXon_tables/" + TaXon_table_xlsx.stem + "_itt.xlsx")
    new_df.to_excel(output_pdf, index=False)

    sg.Popup('All ITT traits have been added to the Taxon table.')


def itt_stats(path_to_outdirs, metadata_to_test, taxonomic_level):

    ## define level to test
    taxonomic_level = 'OTUs'
    meta_data_to_test = 'larva_diet_types_zoophagus'

    ## load TaxonTable
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx).fillna('nan')
    ## collect samples and metadata
    samples = strip_metadata(TaXon_table_df).columns.tolist()[10:]
    categories = set(TaXon_table_df[meta_data_to_test].values.tolist())

    ## adjust taxonomic level if neccessary
    taxon_title = taxonomic_level
    if taxonomic_level in ["ASVs", "ESVs", "OTUs", "zOTUs"]:
        taxon_title = taxonomic_level
        taxonomic_level = "ID"

    ## reduce df to samples, taxon and metadata
    TaXon_table_df = reduce_taxontable(TaXon_table_df, taxonomic_level, samples, meta_data_to_test)

    ## loop through categories and collect number of occurrences per category
    categories_dict = {}
    for category in categories:
        y_sample_counts = []
        for sample in samples:
            present_categories = [i[2] for i in TaXon_table_df[[taxonomic_level, sample, meta_data_to_test]].values.tolist() if i[1] != 0]
            y_sample_counts.append(present_categories.count(category))
        categories_dict[category] = y_sample_counts

    ## loop through samples and collect the number of taxa
    n_taxa_total = []
    for sample in samples:
        n_taxa_total.append(len([i for i in TaXon_table_df[[taxonomic_level, sample]].values.tolist() if i[1] != 0]))

    ## initiate a new figure
    fig = go.Figure()

    ## calculate relative abundance per category
    for key,values in categories_dict.items():
        rel_values = [j/n_taxa_total[i]*100 for i,j in enumerate(values)]

        if key == 'nan':
            key = 'not specified'
        fig = fig.add_trace(go.Bar(x=samples, y=rel_values, name=str(key)))


    fig.update_layout(barmode='stack', width=2000, height=800, template='plotly_white', title=str(meta_data_to_test))
    fig.update_xaxes(tickmode='linear')
    fig.update_yaxes(title=taxon_title + ' (%)')
    fig.show()


#
