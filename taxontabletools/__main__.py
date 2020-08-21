try:
    import glob, sys, os, pkgutil, ast
    import PySimpleGUI as sg
    import matplotlib.pyplot as plt
    from pathlib import Path
    import webbrowser
except:
    print("\nError: You need to install the required packages first!\n")
    exit()

##########################################################################################################################
# update version here (will be displayed on the main layout)
# Support for: u = ubuntu, w = windows, m = macintosh
taxon_tools_version = "Version 1.0.12"

##########################################################################################################################
# general functions

# slice function for lists to split up lists
def slices(list, slice):
    for i in range(0, len(list), slice):
        yield list[i : i + slice]

##########################################################################################################################
##########################################################################################################################

def metadata_table_help_text():
    metadata_table_help_text = """
    The metadata table will be automatically created from the sample names.
    Metadata information is by default separated by an underscore ("_").

    ____________________________________________________________________________

    Example in the newly created metadata table:

     Sample name: ProjectA_No1_RiverX_RegionY

     Raw output name:
            1          2         3          4         (column as metadata)
       ProjectA | No1 | RiverX | RegionY      (rows for samples)

    ____________________________________________________________________________

    Please rename the columns accordingly. This will enhance the recognition
    when choosing a metadata for downstream analysis.
    Each value can be adjusted to the users need before loading the metadata table.

     User defined names:
       Project  |  Sample | Site | Region   (column as metadata)
       ProjectA | No1 | RiverX | RegionY      (rows for samples)

    """
    return metadata_table_help_text

def alpha_diversity_help_text():
    alpha_diversity_help_text = """

    Alpha diversity is calculated with the python package scikit-bio (http://scikit-bio.org/)
    and is based on the number of OTUs per sample, which are displayed as a scatter plot.

    """
    return alpha_diversity_help_text

def beta_diversity_help_text():
    beta_diversity_help_text = """

    Beta diversity is calculated with the python package scikit-bio (http://scikit-bio.org/)
    and is calculated as jaccard-distances, which are illustrated in a distance matrix.
    An R-square test for reliability and validity is conducted.
    Check out matplotlib colormaps to individualize the plot.

    The TaXon table must be converted to incidence data.

    """
    return beta_diversity_help_text

def cca_analysis_help_text():
    cca_analysis_help_text = """

    Canonical-correlation analysis (CCA) analysis plots are calculated with the python package scikit-bio (http://scikit-bio.org/).
    The two axes to plot can be chosen from all available axes.
    The meta data must be a vector of random variables (X = X1...,Xn) and must differ between samples!
    All eigenvalues are saved as an Excel sheet.

    The TaXon table must be converted to incidence data.

    """
    return cca_analysis_help_text

def pcoa_analysis_help_text():
    pcoa_analysis_help_text = """

    Principle coordinate analysis (PCoA) plots are calculated with the python package scikit-bio (http://scikit-bio.org/).
    and are based on jaccard-distances.
    The meta data must differ between samples!
    An R-square test for reliability and validity is conducted.
    The two axes to plot can be chosen from all available axes.
    All eigenvalues are saved as an Excel sheet.

    The TaXon table must be converted to incidence data.

    """
    return pcoa_analysis_help_text

##########################################################################################################################

def main():

    sg.ChangeLookAndFeel('Reddit')

    ##########################################################################################################################
    # start Popup window
    # assign picture paths

    # collect user data from the user_data.txt file in the main TTT directory
    # collect path to __main__.py
    path_to_ttt = os.path.realpath(__file__)
    path_to_ttt = Path(path_to_ttt).parent

    user_data_txt = str(path_to_ttt) + "/user_data.txt"
    f = open(user_data_txt)
    projects_main_path = f.read()

    crash = Path(str(path_to_ttt) + "/_source/crash.png")
    fließgewaesserbewertungde = Path(str(path_to_ttt) + "/_source/fließgewässerbewertung_de.png")
    github = Path(str(path_to_ttt) + "/_source/github.png")
    ttt_main = Path(str(path_to_ttt) + "/_source/main.png")
    ttt_logo = Path(str(path_to_ttt) + "/_source/taxontabletools.png")
    twitter = Path(str(path_to_ttt) + "/_source/twitter.png")

    # fresh start: there is an empty user_data file
    # ask for user Input
    # stay open until a path was defined
    # then write it the user_data file to reload
    while projects_main_path == "":
        projects_main_path = sg.PopupGetFolder("Enter path to ouput directory:", title="Output directory")
        if projects_main_path == None:
            sys.exit()
        f = open(user_data_txt, "w")
        f.write(projects_main_path)
        f.close()

    # create display text
    current_path = "Current path: " + str(projects_main_path)
    # load all available projects
    projects = glob.glob(str(projects_main_path) + '/Projects/*')
    projects_list = []

    for project in projects:
        projects_list.append(Path(project).stem)

    projects_radio = list(slices([sg.Radio(name, "projects", default=True) for name in sorted(projects_list)], 3))

    start_window_layout = [
                [sg.Text('',size=(1,1))],
    			[sg.Text('Output directory', size=(50,1), font=('Arial', 11, "bold"))],
                [sg.Text(current_path)],
                [sg.Text('Define new output directory:')],
                [sg.Input(key = "new_projects_main_path", size=(40,1)), sg.FolderBrowse(), sg.Button("Refresh")],
                [sg.Text('',size=(1,1))],
    			[sg.Text('Project management', size=(50,1), font=('Arial', 11, "bold"))],
                [sg.Text('Create new project folder:')],
                [sg.Input('', key='new_project_folder', size=(40,1)), sg.Button('Create new')],
                [sg.Text('Load existing project folder:')],
                [sg.Frame(layout = projects_radio, title = '')],
                [sg.Button('Load')],
                [sg.Text('',size=(1,1))],
                [sg.Button('Exit', button_color=('black', 'red'))],
                ]

    start_window = sg.Window('Projects', start_window_layout)
    event, values = start_window.read()

    while True:

        new_projects_main_path = values["new_projects_main_path"]

        if event == 'Create new':
            if values["new_project_folder"] != '':
                project_folder = values["new_project_folder"].replace(" ", "_")
                break
            else:
                project_folder = "Default_project"
                break

        if event == 'Load':
            project_folder = ''
            for key, value in values.items():
                if value == True:
                    project_folder = sorted(projects_list)[key]
            if project_folder == '':
                project_folder = "Default_project"
            break

        if event == 'Refresh':
            if new_projects_main_path == None:
                break
            f = open(user_data_txt, "w")
            f.write(new_projects_main_path)
            f.close()
            sg.Popup("Please reload TaxonTableTools to apply changes", title="Refresh output directory")
            sys.exit()

        if event == 'Exit':
            sys.exit()

    start_window.close()

    ##########################################################################################################################
    # check folders and create new folders if neccessary
    try:
        if not os.path.exists(Path(str(projects_main_path) + "/Projects")):
            os.mkdir(Path(str(projects_main_path) + "/Projects"))
    except:
        sg.PopupError("The output directory does not exist anymore! Please refresh the output folder.")
        sys.exit()

    path_to_outdirs = Path(str(projects_main_path) + "/Projects/" + project_folder)
    if not os.path.exists(path_to_outdirs):
        os.mkdir(path_to_outdirs)

    directories_to_create = ["Pie_charts", "Venn_diagrams","TaXon_tables", "TaXon_tables_per_sample", "Taxon_lists", "Rarefaction_curves",
    "Taxonomic_richness_plots", "CCA_plots", "Taxonomic_resolution_plots", "Meta_data_table", "Site_occupancy_plots", "Read_proportions_plots", "Basic_stats",
    "Krona_charts", "Perlodes", "Alpha_diversity", "Beta_diversity", "PCoA_plots"]

    for directory in directories_to_create:
        dirName = Path(str(path_to_outdirs) + "/" + directory + "/")
        if not os.path.exists(dirName):
            os.mkdir(Path(dirName))

    ##########################################################################################################################
    ##########################################################################################################################
    ##########################################################################################################################

    data_conversion = [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('TaXon table converter', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('Enter path to taxonomy table and read table to create a new table in TaXon format', size=(85,2))],
    					[sg.Text('Taxonomy table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'taxonomy_table_path'), sg.Button("Check table", key = 'run_taxonomy_table_check', button_color=('black', 'white'))],
                        [sg.Text('Sorting method:', size=(20, 1)),  sg.Radio('First hit', "sorting_method", key='sorting_method_fh'), sg.Radio('JAMP', "sorting_method", default=True, key='sorting_method_jamp'), sg.Radio('BOLDigger', "sorting_method", key='sorting_method_boldigger')],
    					[sg.Text('Read table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'Read_table_path'), sg.Button("Check table", key = 'run_read_table_check', button_color=('black', 'white'))],
                        [sg.Text('Format:', size=(20, 1)), sg.Radio('TTT', "read_table_format", key='read_table_format_TTT', default=True),  sg.Radio('JAMP', "read_table_format", key='read_table_format_jamp'), sg.Radio('QIIME2', "read_table_format", key='read_table_format_qiime')],
                        [sg.Text('Name of TaXon table:', size=(20, 1)), sg.Input(project_folder + "_taxon_table", key = 'TaXon_table_name'), sg.Text('.xlsx'), sg.Button("Convert", key = 'run_taxon_converter')],
    					[sg.Text('')],
                        [sg.Text('Check TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'taXon_table_to_check'), sg.Button("Check table", key = 'run_taXon_table_check', button_color=('black', 'white'))],
                        [sg.Text('_'*105)],
    					[sg.Image(ttt_main)],
    					]

    processing_1_layout = [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Processing #1', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'taxon_table_processing_1_table')],
    					[sg.Text('_'*105)],
    					[sg.Text('Replicates',size=(40,2), font=('Arial', 11, "bold"))],
                        [sg.Text("Replicate suffixes:", size=(30, 1)), sg.Input("a, b", size=(15,1), key="replicate_suffix_entry"),
                        sg.Text("", size=(1,1)), sg.Button("Preview", key="run_preview_suffix"), sg.Text("", size=(1,1)),
                        sg.Frame(layout=[[sg.Text("Comma separated list of suffixes.")]], title="Information")],
                        [sg.Text('',size=(1,1))],
    					[sg.Text("Replicate consistency", size=(30, 1)), sg.Button("Run", key='run_replicate_consistency_filter'), sg.Text("", size=(1,1)),
                        sg.Frame(layout=[[sg.Text("Keep OTUs present in all replicates.")]], title="Information")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text("Replicate merging", size=(30, 1)), sg.Button("Run", key='run_combine_replicates'), sg.Text("", size=(1,1)),
                        sg.Frame(layout=[[sg.Text("Merge replicates to one entry.")]], title="Information")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('TaXon table per sample', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Create TaXon tables per-sample", size=(30,1)), sg.Button("Run", key = 'run_taxon_table_per_sample')],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Presence absence table', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Convert to presence absence TaXon table", size=(30,1)), sg.Button("Run", key = 'run_convert_to_pa')],
    					[sg.Text('',size=(1,1))],
    					]

    processing_2_layout = [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Processing #2', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'taxon_table_processing_2_table')],
    					[sg.Text('_'*105)],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Taxon-based filter',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Taxonomic level: "), sg.Radio('Phylum', "tt_filter", default=True, key='filter_phylum'), sg.Radio('Class', "tt_filter", key='filter_class'),
    					sg.Radio('Order', "tt_filter", key='filter_order'), sg.Radio('Family', "tt_filter", key='filter_family'),
    					sg.Radio('Genus', "tt_filter", key='filter_genus'), sg.Radio('Species', "tt_filter", key='filter_species')],
    					[sg.Text("", size=(20, 1))],
    					[sg.Text("Taxon-based filter:", size=(25, 1)), sg.Button("Run", key = 'run_taxon_filtering'), sg.Text("", size=(2,1)), sg.Frame(layout=[[sg.Radio('all', "cm_taxon_based_filter", key='cm_taxon_based_filter_all', default=True),
                        sg.Radio('none', "cm_taxon_based_filter", key='cm_taxon_based_filter_none')]], title="Check marks")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Sample-based filter', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Sample-based filter:", size=(25, 1)), sg.Button("Run", key = 'run_sample_filtering'), sg.Text("", size=(2,1)), sg.Frame(layout=[[sg.Radio('all', "cm_sample_based_filter", key='cm_sample_based_filter_all', default=True),
                        sg.Radio('none', "cm_sample_based_filter", key='cm_sample_based_filter_none')]], title="Check marks")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Metadata table', size=(50,2), font=('Arial', 11, "bold")), sg.Text("", size=(2,1))],
    					[sg.Text("Create metadata table", size=(25,1)), sg.Button("Run", key = 'run_create_meta_data_table'), sg.Text("", size=(1,1)),
                        sg.Button("How does it work?", button_color=('black', 'white'), key="run_metadata_table_help" )],
    					[sg.Text('',size=(1,1))],
    					]

    analysis_1_layout = [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('TaXon table analyses #1', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'statistics_1_taxon_table_path')],
    					[sg.Text('_'*105)],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Basic statistics',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Calculate basic statistics", size=(25,1)), sg.Button("Run", key = 'run_basic_stats'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("8", size=(3,1), key="bstats_w"),
                        sg.Input("8", size=(3,1), key="bstats_h"), sg.Text("", size=(1,1)), sg.Text("Font size:"), sg.Input("5", size=(3,1), key="bstats_f")]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Taxonomic resolution',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Taxonomic resolution plot", size=(25,1)), sg.Button("Run", key = 'run_taxonomic_resolution'), sg.Text("", size=(1,1)),
                        sg.Frame(layout=[
                        [sg.Text("Plot size (w,h, ylim):"), sg.Input("8", size=(3,1), key="x_tax_res"), sg.Input("8", size=(3,1), key="y_tax_res"), sg.Input("", size=(4,1), key="ylim_tax_res"), sg.Text("", size=(1,1)), sg.Text("Font size:"), sg.Input("10", size=(3,1), key="font_tax_res")],
                        [sg.Text("Plot type:"), sg.Radio('a', "tres_plot_type", key='tres_type_a', default=True), sg.Radio('b', "tres_plot_type", key='tres_type_b')]
                        ], title="Settings", size=(1,2))],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Taxonomic richness',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Taxonomic richness plot", size=(25,1)), sg.Button("Run", key = 'run_taxonomic_richness'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h,ylim):"), sg.Input("8", size=(3,1), key="x_tax_rich"),
                        sg.Input("8", size=(3,1), key="y_tax_rich"), sg.Input("", size=(4,1), key="ylim_tax_rich"), sg.Text("", size=(1,1)), sg.Text("Font size:"), sg.Input("10", size=(3,1), key="font_tax_rich")]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					]

    analysis_2_layout =  [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('TaXon table analyses #2', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'statistics_2_taxon_table_path')],
    					[sg.Text('_'*105)],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('OTU abundance pie charts',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Create pie charts", size=(20,1)), sg.Button("Run", key = 'run_pie_chart'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Text("8", size=(3,1)),
                        sg.Text("8", size=(3,1)), sg.Text("", size=(1,1)), sg.Text("Font size:"), sg.Input("5", size=(3,1), key="pc_label_font_size")]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Venn diagrams',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text('Compare to table(s):', size=(16, 1)), sg.Input("", size=(10,1)), sg.FileBrowse(key = 'venn_taxon_table_2_path'), sg.Input("", size=(10,1)), sg.FileBrowse(key = 'venn_taxon_table_3_path'),
                        sg.Text("", size=(1,1)), sg.Text("Output file:"), sg.Input("venn_diagram", size=(15,1), key="venn_diagram_name")],
    					[sg.Text("", size=(1,1))],
    					[sg.Text("Compare TaXon tables", size=(20,1)), sg.Button("venn2", key = 'run_venn2_diagram'), sg.Button("venn3", key = 'run_venn3_diagram'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Text("8", size=(3,1)),
                        sg.Text("8", size=(3,1)), sg.Text("", size=(1,1)), sg.Text("Font size:"), sg.Text("10", size=(3,1))]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Rarefaction curve',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Calculate rarefaction curve with"), sg.Input('1000', size=(6,1), key='repetitions'), sg.Text('repetitions'),
                        sg.Frame(layout=[[sg.Text("Y-limit:"), sg.Input("", size=(4,1), key="rarefaction_ylim"), sg.Text("Error bar style:"), sg.Radio('a', "error_bar_style", key='error_bar_style_a', default=True), sg.Radio('b', "error_bar_style", key='error_bar_style_b')]], title="Settings")],
                        [sg.Text("OTUs: "), sg.Button("Run", key = 'run_rarefaction_curve_OTUs'), sg.Text("", size=(1,1)), sg.Text("Species: "), sg.Button("Run", key = 'run_rarefaction_curve_species'), sg.Text("", size=(1,1))],
    					[sg.Text('',size=(1,1))],
    					]

    analysis_3_layout =  [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('TaXon table analyses #3', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'statistics_3_taxon_table_path')],
    					[sg.Text('_'*105)],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Read proportions',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Taxonomic level: "), sg.Radio('OTUs', "ra_filter", key='ra_OTUs'), sg.Radio('Phylum', "ra_filter", key='ra_phylum'), sg.Radio('Class', "ra_filter", key='ra_class'),
    					sg.Radio('Order', "ra_filter", key='ra_order'), sg.Radio('Family', "ra_filter", key='ra_family'),
    					sg.Radio('Genus', "ra_filter", key='ra_genus'), sg.Radio('Species', "ra_filter", default=True, key='ra_species')],
    					[sg.Text("", size=(20, 1))],
    					[sg.Text("Calculate read proportions", size=(25,1)), sg.Button("Run", key = 'run_read_proportions'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h): "), sg.Input("8", size=(3,1), key="x_read_props_scatter"),
                        sg.Input("8", size=(3,1), key="y_read_props_scatter"), sg.Text("", size=(1,1)), sg.Text("Font size:"), sg.Input("4", size=(3,1), key="font_read_props_scatter")],
                        [sg.Text("No hit visualization:"), sg.Radio('as nan', "read_props_taxon_level", key='read_props_as_nan', default=True), sg.Radio('best hit', "read_props_taxon_level", key='read_props_best_hit')],
                        [sg.Checkbox("Alternating colors", key="read_props_alternating_colors", default=True)]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Site occupancy', size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Calculate site occupancy", size=(25,1)), sg.Button("Run", key = 'run_site_occupancy'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h): "), sg.Input("8", size=(3,1), key="x_site_occ"),
                        sg.Input("8", size=(3,1), key="y_site_occ"), sg.Text("", size=(1,1)), sg.Text("Font size:"), sg.Input("8", size=(3,1), key="font_site_occ")]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Krona charts',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Create Krona chart", size=(25,1)), sg.Button("Run", key = 'run_create_krona_chart'), sg.Text("", size=(1,1)),
                        sg.Frame(layout=[[sg.Text("Charts are based on read proportions. Krona tools is required.")]], title="Information")],
    					[sg.Text('',size=(1,1))],
    					]

    analysis_4_layout =  [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('TaXon table analyses #4', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'statistics_4_taxon_table_path')],
    					[sg.Text('_'*105)],
    					[sg.Text('',size=(1,1))],

    					[sg.Text('Diversity analyses',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Calculate alpha diversity", size=(23,1)), sg.Button("Run", key = 'run_alpha_diversity'),
                        sg.Button("Help", key = 'run_alpha_diversity_help_text', button_color=('black', 'white')),
                        sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("8", size=(3,1), key="alpha_w"),
                        sg.Input("8", size=(3,1), key="alpha_h"), sg.Text("Font size:"), sg.Input("8", size=(3,1), key="alpha_font")],
                        [sg.Radio("Scatter plot", "alpha_plot", key="alpha_scatter_plot", default=True), sg.Text("Scatter size:"), sg.Input("40", size=(3,1), key="alpha_s")],
                        [sg.Radio("Box plot", "alpha_plot", key="alpha_box_plot")]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text("Calculate beta diversity", size=(23,1)), sg.Button("Run", key = 'run_beta_diversity'),
                        sg.Button("Help", key = 'run_beta_diversity_help_text', button_color=('black', 'white')),
                        sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("8", size=(3,1), key="beta_w"),
                        sg.Input("8", size=(3,1), key="beta_h"), sg.Text("Font size:"), sg.Input("8", size=(3,1), key="beta_font")],
                        [sg.Text("Cmap:"), sg.Input("Blues_r", size=(15,1), key="beta_cmap")]], title="Settings")],
    					[sg.Text('',size=(1,1))],

    					[sg.Text('Ordination analyses',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Calculate CCA", size=(23,1)), sg.Button("Run", key = 'run_CCA_analysis'),
                        sg.Button("Help", key = 'run_cca_analysis_help_text', button_color=('black', 'white')),
                        sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("8", size=(3,1), key="cca_w"),
                        sg.Input("8", size=(3,1), key="cca_h"), sg.Text("", size=(1,1)), sg.Text("Scatter size:"), sg.Input("10", size=(3,1), key="cca_s")]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text("Calculate PCoA", size=(23,1)), sg.Button("Run", key = 'run_PCoA_analysis'),
                        sg.Button("Help", key = 'run_pcoa_analysis_help_text', button_color=('black', 'white')),
                        sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("8", size=(3,1), key="pcoa_w"),
                        sg.Input("8", size=(3,1), key="pcoa_h"), sg.Text("", size=(1,1)), sg.Text("Scatter size:"), sg.Input("10", size=(3,1), key="pcoa_s")]], title="Settings")],
    					[sg.Text('',size=(1,1))],
        				]

    create_taxon_list_layout = [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Create taxon list',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'taxon_list_taxon_table_path')],
    					[sg.Text('_'*105)],
                        [sg.Text("Enter name of output file:"), sg.Input(project_folder + "_taxon_list", size=(48,1), key='taxon_list_output_file_name'), sg.Text(".txt + .xlsx")],
                        [sg.Text('Choose output language:'), sg.Radio('English', 'Language', key='English', default=True), sg.Radio('German', 'Language', key='German', default=False)],
                        [sg.Text('GBIF link:', size=(20,1)), sg.CB("", key="create_gbif_link", default=True), sg.Frame(layout=[[sg.Text("Internet connection is required.")]], title="Information")],
                        [sg.Text('Intraspecific distances:', size=(20,1)), sg.CB("", key="calc_dist", default=True), sg.Frame(layout=[[sg.Text("This option can drastically increase run times.")]], title="Information")],
                        [sg.Text('Occupancy per sample:', size=(20,1)), sg.CB("", key="calc_occupancy", default=True), sg.Frame(layout=[[sg.Text("This option can drastically increase run times.")]], title="Information")],
                        [sg.Text('', size=(1,1))],
                        [sg.Text('Additional information', font=('Arial', 10, "bold"))],
                        [sg.Text('Description:', size=(20,1)), sg.Input('', size=(48,1), key='TL_description')], # TL = taxon list -> makes it easier to find
                        [sg.Text('Author(s):', size=(20,1)), sg.Input('', size=(48,1), key='TL_authors')],
                        [sg.Text('Lab protocol:', size=(20,1)), sg.Input('', size=(48,2), key='TL_lab_protocol')],
                        [sg.Text('No. of replicates:', size=(20,1)), sg.Input('X (extraction); X (1-st-step); X (2nd-step)', size=(48,1), key='TL_replicates')],
                        [sg.Text('No. of neg. contr.:', size=(20,1)), sg.Input('', size=(48,1), key='TL_negative_controls')],
                        [sg.Text('Primers:', size=(20,1)), sg.Input('', size=(48,1), key='TL_primers')],
                        [sg.Text('Sequencing run:', size=(20,1)), sg.Input('', size=(48,1), key='TL_sequencing_run')],
                        [sg.Text('Bionf. pipeline:', size=(20,1)), sg.Input('', size=(48,1), key='TL_pipeline')],
                        [sg.Text('', size=(1,1))],
                        [sg.Text("Create taxon list"), sg.Button("Run", key = 'run_create_taxon_list')],
    					]

    create_perlodes_input_layout = [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Water Framework Directive (WFD)',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'perlodes_taxon_table_path')],
    					[sg.Text('_'*105)],
    					[sg.Text('',size=(1,1))],
        				[sg.Text('Create Perlodes input file',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text('Operational taxon list:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'operational_taxon_list_path'),
                        sg.Button('Download', key='open_operational_taxonlist_download', button_color=('black', 'white'))],
    					[sg.Text('',size=(1,1))],
    					[sg.Text("Create Perlodes input file: "), sg.Button("Run", key = 'run_create_perlodes_input_file')],
    					[sg.Text('',size=(1,2))],
    					[sg.Button(key = 'open_fgbewertung', button_color=('white', 'white'), image_filename=fließgewaesserbewertungde, image_size=(687,157))],
                        ]

    layout = [  [sg.Image(ttt_logo), sg.Text("", size=(9,1)), sg.Text('Project:', font=('Arial', 12, "bold")), sg.Text(project_folder, font=('Arial', 12, "bold"))],
    			[sg.Text('',size=(1,1))],
    			[sg.TabGroup([[
                sg.Tab('Data conversion', data_conversion),
                sg.Tab('Processing #1', processing_1_layout),
                sg.Tab('Processing #2', processing_2_layout),
                sg.Tab('Analysis #1', analysis_1_layout),
                sg.Tab('Analysis #2', analysis_2_layout),
                sg.Tab('Analysis #3', analysis_3_layout),
                sg.Tab('Analysis #4', analysis_4_layout),
                sg.Tab('Taxon list', create_taxon_list_layout),
                sg.Tab('WFD', create_perlodes_input_layout)]])],
    			[sg.Text('',size=(1,1))],
    			[sg.Exit(button_color=('black', 'red')), sg.Text("", size=(57,1)),
                sg.Button("Log file network", button_color=('black', 'white'), key='run_log_network'), sg.Text("", size=(4,1)),
                sg.Button(image_filename=github, image_size=(26,26), key='open_github', button_color=('black', 'white')),
                sg.Button(key='open_twitter', button_color=('white', 'white'), image_filename=twitter, image_size=(26,26)),
                sg.Text('', size=(1,1)), sg.Text(taxon_tools_version, font=('Arial', 8))]]


    # Create the Window
    window = sg.Window('TaxonTableTools', layout, grab_anywhere=True)
    win2_active=False

    ##########################################################################################################################

    while True:
        try:

            event, values = window.Read()

            if event is None or event == 'Exit':
                break

            i = 0
            # define variables
            TaXon_table_name = values['TaXon_table_name']
            taxon_table_processing_1_table = values['taxon_table_processing_1_table']
            replicate_suffix_entry = values['replicate_suffix_entry']
            taxonomy_table_path = values['taxonomy_table_path']
            Read_table_path = values['Read_table_path']
            venn_taxon_table_2_path = values['venn_taxon_table_2_path']
            venn_taxon_table_3_path = values['venn_taxon_table_3_path']
            venn_diagram_name = values['venn_diagram_name']
            taxon_list_taxon_table_path = values['taxon_list_taxon_table_path']
            taxon_table_processing_2_table = values['taxon_table_processing_2_table']
            filter_phylum = values['filter_phylum']
            filter_order = values['filter_order']
            filter_genus = values['filter_genus']
            filter_class = values['filter_class']
            filter_family = values['filter_family']
            filter_species = values['filter_species']
            ra_OTUs = values['ra_OTUs']
            ra_phylum = values['ra_phylum']
            ra_order = values['ra_order']
            ra_genus = values['ra_genus']
            ra_class = values['ra_class']
            ra_family = values['ra_family']
            ra_species = values['ra_species']
            statistics_1_taxon_table_path = values['statistics_1_taxon_table_path']
            statistics_2_taxon_table_path = values['statistics_2_taxon_table_path']
            statistics_3_taxon_table_path = values['statistics_3_taxon_table_path']
            statistics_4_taxon_table_path = values['statistics_4_taxon_table_path']
            taxon_list_output_file_name = values['taxon_list_output_file_name']
            English = values['English']
            German = values['German']
            repetitions = values['repetitions']
            x_tax_res = values["x_tax_res"]
            y_tax_res = values["y_tax_res"]
            ylim_tax_res = values["ylim_tax_res"]
            font_tax_res = values["font_tax_res"]
            x_tax_rich = values["x_tax_rich"]
            y_tax_rich = values["y_tax_rich"]
            ylim_tax_rich = values["ylim_tax_rich"]
            font_tax_rich = values["font_tax_rich"]
            x_read_props_scatter = values["x_read_props_scatter"]
            y_read_props_scatter = values["y_read_props_scatter"]
            font_read_props_scatter = values["font_read_props_scatter"]
            x_site_occ = values["x_site_occ"]
            y_site_occ = values["y_site_occ"]
            font_site_occ = values["font_site_occ"]
            error_bar_style_a = values["error_bar_style_a"]
            error_bar_style_b = values["error_bar_style_b"]
            operational_taxon_list_path = values["operational_taxon_list_path"]
            perlodes_taxon_table_path = values["perlodes_taxon_table_path"]
            alpha_w = values["alpha_w"]
            alpha_h = values["alpha_h"]
            alpha_s = values["alpha_s"]
            alpha_font = values["alpha_font"]
            beta_w = values["beta_w"]
            beta_h = values["beta_h"]
            beta_font = values["beta_font"]
            beta_cmap = values["beta_cmap"]
            tres_type_a = values["tres_type_a"]
            tres_type_b = values["tres_type_b"]
            cm_taxon_based_filter_all = values['cm_taxon_based_filter_all']
            cm_taxon_based_filter_none = values['cm_taxon_based_filter_none']
            cm_sample_based_filter_all = values['cm_sample_based_filter_all']
            cm_sample_based_filter_none = values['cm_sample_based_filter_none']
            create_gbif_link = values['create_gbif_link']
            calc_dist = values['calc_dist']
            calc_occupancy = values['calc_occupancy']
            pc_label_font_size = values['pc_label_font_size']
            sorting_method_fh = values['sorting_method_fh']
            sorting_method_jamp = values['sorting_method_jamp']
            sorting_method_boldigger = values['sorting_method_boldigger']
            rarefaction_ylim = values['rarefaction_ylim']
            cca_w = values['cca_w']
            cca_h = values['cca_h']
            cca_s = values['cca_s']
            pcoa_w = values['pcoa_w']
            pcoa_h = values['pcoa_h']
            pcoa_s = values['pcoa_s']
            bstats_h = values["bstats_h"]
            bstats_w = values["bstats_w"]
            bstats_f = values["bstats_f"]
            read_props_as_nan = values["read_props_as_nan"]
            read_props_best_hit = values["read_props_best_hit"]
            taXon_table_to_check = values["taXon_table_to_check"]
            read_table_format_TTT = values["read_table_format_TTT"]
            read_table_format_jamp = values["read_table_format_jamp"]
            read_table_format_qiime = values["read_table_format_qiime"]
            read_props_alternating_colors = values["read_props_alternating_colors"]
            alpha_box_plot = values["alpha_box_plot"]
            alpha_scatter_plot = values["alpha_scatter_plot"]

            print("########", "\n")

            if event == 'run_taxon_converter':

                if (taxonomy_table_path == '' or Read_table_path == ''):
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    from taxontabletools.taxon_converter import taxon_converter
                    if sorting_method_fh == True:
                        sheet_name = "First hit"
                    elif sorting_method_jamp == True:
                        sheet_name = "JAMP hit"
                    else:
                        sheet_name = "BOLDigger hit"
                    taxon_converter(Read_table_path, taxonomy_table_path, TaXon_table_name, sheet_name, path_to_outdirs)

            if event == 'run_taxonomy_table_check':
                print("Checking taxonomy table format")

                if (taxonomy_table_path == ''):
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    from taxontabletools.check_taxononomy_table_format import check_taxononomy_table_format
                    if sorting_method_fh == True:
                        sheet_name = "First hit"
                    elif sorting_method_jamp == True:
                        sheet_name = "JAMP hit"
                    else:
                        sheet_name = "BOLDigger hit"
                    check_taxononomy_table_format(taxonomy_table_path, sheet_name)

            if event == 'run_read_table_check':
                print("Checking read table format")

                if (Read_table_path == ''):
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                elif read_table_format_TTT == True:
                    from taxontabletools.check_read_table_format import check_read_table_format_TTT
                    check_read_table_format_TTT(Read_table_path)
                elif read_table_format_jamp == True:
                    sg.Popup("Not supported yet!")
                elif read_table_format_qiime == True:
                    sg.Popup("Not supported yet!")

            if event == 'run_taXon_table_check':
                print("Checking TaXon table format")

                if (taXon_table_to_check == ''):
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    from taxontabletools.check_taXon_table_format import check_taXon_table_format
                    check_taXon_table_format(taXon_table_to_check)

            if event == "run_taxon_table_per_sample":
                print("Creating OTU table for each sample")

                if taxon_table_processing_1_table == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    from taxontabletools.create_taxon_table_per_sample import create_taxon_table_per_sample
                    create_taxon_table_per_sample(taxon_table_processing_1_table, path_to_outdirs)

            if event == "run_convert_to_pa":
                print("Converting to presence absence table")

                if taxon_table_processing_1_table == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    from taxontabletools.convert_to_presence_absence import convert_to_presence_absence
                    convert_to_presence_absence(taxon_table_processing_1_table, path_to_outdirs)

            if event == 'run_pie_chart':
                print("Pie chart")

                if statistics_2_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    from taxontabletools.pie_charts import pie_charts
                    pie_charts(statistics_2_taxon_table_path, path_to_outdirs, pc_label_font_size)

            if event == 'run_venn2_diagram':
                print("Venn diagrams")

                if os.path.exists(str(path_to_outdirs) + "/Venn_diagrams/" + venn_diagram_name):
                    sg.PopupError("Error: Venn directory already exists", keep_on_top=True)
                    print("Error: Venn directory already exists")
                elif (statistics_2_taxon_table_path == '' or venn_taxon_table_2_path == ''):
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    from taxontabletools.venn_diagram import venn_diagram
                    venn_taxon_table_3_path = False
                    venn_diagram(statistics_2_taxon_table_path, venn_taxon_table_2_path, venn_taxon_table_3_path, venn_diagram_name, path_to_outdirs)

            if event == 'run_venn3_diagram':
                print("Venn diagrams")

                if os.path.exists(str(path_to_outdirs) + "/Venn_diagrams/" + venn_diagram_name):
                    sg.PopupError("Error: Venn directory already exists", keep_on_top=True)
                    print("Error: Venn directory already exists")
                elif (statistics_2_taxon_table_path == '' or venn_taxon_table_2_path == ''):
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    from taxontabletools.venn_diagram import venn_diagram
                    venn_diagram(statistics_2_taxon_table_path, venn_taxon_table_2_path, venn_taxon_table_3_path, venn_diagram_name, path_to_outdirs)

            if event == 'run_rarefaction_curve_OTUs':
                if statistics_2_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    if error_bar_style_a == True:
                        error_bar_style = "a"
                    if error_bar_style_b == True:
                        error_bar_style = "b"
                    from taxontabletools.rarefaction_curve import rarefaction_curve_OTUs
                    rarefaction_curve_OTUs(statistics_2_taxon_table_path, int(repetitions), error_bar_style, rarefaction_ylim, path_to_outdirs)

            if event == 'run_rarefaction_curve_species':
                if statistics_2_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    if error_bar_style_a == True:
                        error_bar_style = "a"
                    if error_bar_style_b == True:
                        error_bar_style = "b"
                    from taxontabletools.rarefaction_curve import rarefaction_curve_species
                    rarefaction_curve_species(statistics_2_taxon_table_path, int(repetitions), error_bar_style, rarefaction_ylim, path_to_outdirs)

            if event == 'run_create_taxon_list':
                if English == True:
                    language = "English"
                if German == True:
                    language = "German"
                if taxon_list_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    from taxontabletools.create_taxon_list import create_taxon_list
                    create_taxon_list(taxon_list_taxon_table_path, taxon_list_output_file_name, language, values, create_gbif_link, calc_dist, calc_occupancy, taxon_tools_version, path_to_outdirs)

            if event == 'run_taxon_filtering' and not win2_active:
                if filter_phylum == True:
                    mask = "Phylum"
                elif filter_class == True:
                    mask = "Class"
                elif filter_order == True:
                    mask = "Order"
                elif filter_family == True:
                    mask = "Family"
                elif filter_genus == True:
                    mask = "Genus"
                elif filter_species == True:
                    mask = "Species"

                if taxon_table_processing_2_table == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")

                else:
                    from taxontabletools.taxon_table_filtering import create_taxa_mask
                    available_taxa = create_taxa_mask(taxon_table_processing_2_table, mask)
                    win2_active = True
                    window.Hide()
                    if cm_taxon_based_filter_all == True:
                        available_taxa_list = list(slices([sg.CB(name, default=True) for name in sorted(available_taxa)], 8))
                    else:
                        available_taxa_list = list(slices([sg.CB(name) for name in sorted(available_taxa)], 8))
                    layout2 = [[sg.Text("Taxon filter", size=(20,1))],
                    [sg.Text('Enter appendix for output file: '), sg.InputText('file_appendix', key = 'appendix_name', size=(35,1))],
                    [sg.Text("Filter taxa below a similarity of: "), sg.Input('0', size=(4,1), key='threshold'), sg.Text("%")],
                    [sg.Button('Filter')],
                    [sg.Button('Back')],
                    [sg.Text("")],
                    [sg.Frame(layout = available_taxa_list, title = 'Check taxa to exlude')],
                    [sg.Text("")]]

                    win2 = sg.Window('TaXon table filtering', layout2, keep_on_top=True, grab_anywhere=True)
                    while True:
                        event2, values2 = win2.Read()

                        if event2 == 'Filter':

                            i = 0
                            check_dict = {}
                            check_list = []
                            appendix_name = values2['appendix_name']
                            threshold = values2['threshold']
                            del values2['appendix_name']
                            del values2['threshold']
                            for key, value in values2.items():
                                check_dict[available_taxa[i]] = value
                                i += 1
                            for key, value in check_dict.items():
                                if value == True:
                                    check_list.append(key)

                            from taxontabletools.taxon_table_filtering import taxon_table_filter

                            if check_list == []:
                                check_list = ['no taxa to filter']
                                taxon_table_filter(taxon_table_processing_2_table, check_list, mask, appendix_name, threshold, path_to_outdirs)

                            else:
                                taxon_table_filter(taxon_table_processing_2_table, check_list, mask, appendix_name, threshold, path_to_outdirs)

                            win2.Close()
                            win2_active = False
                            window.UnHide()
                            break

                        if event2 is None or event2 == 'Back':
                            win2.Close()
                            win2_active = False
                            window.UnHide()
                            break

            if event == 'run_sample_filtering' and not win2_active:
                if taxon_table_processing_2_table == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    from taxontabletools.taxon_table_filtering import create_sample_mask
                    available_samples = sorted(create_sample_mask(taxon_table_processing_2_table))
                    win2_active = True
                    window.Hide()
                    if cm_sample_based_filter_all == True:
                        available_samples_list = list(slices([sg.CB(name, default=True) for name in sorted(available_samples)], 6))
                    else:
                        available_samples_list = list(slices([sg.CB(name) for name in sorted(available_samples)], 6))
                    layout2 = [[sg.Text("Sample filter", size=(20,1))],
                    [sg.Frame(layout = available_samples_list, title = 'Check samples to exlude')],
                    [sg.Text('Enter appendix for output file: '), sg.InputText('file_appendix', key = 'appendix_name', size=(35,1))],
                    [sg.Button('Filter')],
                    [sg.Button('Back')]]

                    win2 = sg.Window('TaXon table filtering', layout2, keep_on_top=True, grab_anywhere=True)

                    while True:
                        event2, values2 = win2.Read()

                        if event2 == 'Filter':
                            i = 0
                            check_dict = {}
                            check_list = []
                            appendix_name = values2['appendix_name']
                            del values2['appendix_name']
                            for key, value in values2.items():
                                check_dict[available_samples[i]] = value
                                i += 1
                            for key, value in check_dict.items():
                                if value == True:
                                    check_list.append(key)

                            print("Excluded samples:", ', '.join(check_list))
                            from taxontabletools.taxon_table_filtering import exclude_samples
                            exclude_samples(taxon_table_processing_2_table, check_list, appendix_name, path_to_outdirs)

                            win2.Close()
                            win2_active = False
                            window.UnHide()
                            break

                        if event2 is None or event2 == 'Back':
                            win2.Close()
                            win2_active = False
                            window.UnHide()
                            break

            if event == 'run_combine_replicates':

                replicate_suffix_list = replicate_suffix_entry.replace(" ", "").split(",")

                if taxon_table_processing_1_table == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    from taxontabletools.combine_replicates import combine_replicates
                    combine_replicates(taxon_table_processing_1_table, replicate_suffix_list, path_to_outdirs)

            if event == 'run_preview_suffix':
                replicate_suffix_list = replicate_suffix_entry.replace(" ", "").split(",")
                text = "\n" + "Sample_No1_comb = " + ', '.join(["Sample_No1_" + suffix for suffix in replicate_suffix_list]) + "\n"
                sg.Popup(text, title='Suffix preview')

            if event == 'run_taxonomic_richness':
                if statistics_1_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    from taxontabletools.calculate_taxonomic_richness import calculate_taxonomic_richness
                    calculate_taxonomic_richness(statistics_1_taxon_table_path, path_to_outdirs, x_tax_rich, y_tax_rich, font_tax_rich, ylim_tax_rich)

            if event == 'run_taxonomic_resolution':
                if statistics_1_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    if tres_type_a == True:
                        figure_type = "a"
                    if tres_type_b == True:
                        figure_type = "b"
                    from taxontabletools.calculate_taxonomic_resolution import calculate_taxonomic_resolution
                    calculate_taxonomic_resolution(statistics_1_taxon_table_path, path_to_outdirs, x_tax_res, y_tax_res, font_tax_res, figure_type, ylim_tax_res)

            if event == 'run_create_meta_data_table':
                if taxon_table_processing_2_table == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    from taxontabletools.create_metadata_table import create_metadata_table
                    create_metadata_table(taxon_table_processing_2_table, path_to_outdirs)

            if event == 'run_site_occupancy':
                if statistics_3_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")

                else:
                    from taxontabletools.get_available_meta_data import get_available_meta_data
                    meta_data_to_test = get_available_meta_data(statistics_3_taxon_table_path, path_to_outdirs)

                    if meta_data_to_test == False:
                                 sg.PopupError("Metadata file not found!", keep_on_top=True)
                                 print("Metadata file not found!")
                                 print("\n" + "########")
                    else:
                                 win2_active = True
                                 window.Hide()
                                 meta_data_to_test_list = list(slices([sg.Radio(name, "metadata", key=name, default=True) for name in sorted(meta_data_to_test)], 5))
                                 taxonmic_levels_to_choose = ['OTUs', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
                                 taxonmic_levels_list = list(slices([sg.Radio(name, "taxa", key=name, default=True) for name in taxonmic_levels_to_choose], 6))
                                 layout2 = [[sg.Text("Site occupancy", size=(20,1))],
                                 [sg.Frame(layout = meta_data_to_test_list, title = 'Check metadata to test')],
                                 [sg.Text('',size=(1,1))],
                                 [sg.Frame(layout = taxonmic_levels_list, title = 'Check taxonomic level to use')],
                                 [sg.Text('',size=(1,1))],
                                 [sg.Button('Calculate')],
                                 [sg.Button('Back')]] # test = True

                                 win2 = sg.Window('Site occpancy', layout2, keep_on_top=True)

                                 while True:
                                              event2, values2 = win2.Read()

                                              if event2 == 'Calculate':
                                                           for input_value, test in values2.items():
                                                               if (test == True and input_value in meta_data_to_test):
                                                                   meta_data_to_test = input_value
                                                               if (test == True and input_value in taxonmic_levels_to_choose):
                                                                   taxon_level = input_value
                                                           from taxontabletools.site_occupancy import site_occupancy
                                                           site_occupancy(statistics_3_taxon_table_path, meta_data_to_test, taxon_level, path_to_outdirs, x_site_occ, y_site_occ, font_site_occ)
                                                           win2.Close()
                                                           win2_active = False
                                                           window.UnHide()

                                              if event2 is None or event2 == 'Back':
                                                           win2.Close()
                                                           win2_active = False
                                                           window.UnHide()
                                                           break

                                 print("\n" + "########")

            if event == 'run_read_proportions':
                if statistics_3_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    if ra_OTUs == True:
                        taxonomic_level = "OTUs"
                    if ra_phylum == True:
                        taxonomic_level = "Phylum"
                    elif ra_class == True:
                        taxonomic_level = "Class"
                    elif ra_order == True:
                        taxonomic_level = "Order"
                    elif ra_family == True:
                        taxonomic_level = "Family"
                    elif ra_genus == True:
                        taxonomic_level = "Genus"
                    elif ra_species == True:
                        taxonomic_level = "Species"

                    if read_props_as_nan == True:
                        no_hit_visualization = "as_nan"
                    if read_props_best_hit == True:
                        no_hit_visualization = "best_hit"

                    from taxontabletools.read_proportions import read_proportions
                    read_proportions(statistics_3_taxon_table_path, taxonomic_level, path_to_outdirs, x_read_props_scatter, y_read_props_scatter, font_read_props_scatter, no_hit_visualization, read_props_alternating_colors)

            if event == 'run_basic_stats':
                if statistics_1_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    from taxontabletools.basic_stats import basic_stats
                    basic_stats(statistics_1_taxon_table_path, bstats_h, bstats_w, bstats_f, path_to_outdirs)

            if event == 'run_create_krona_chart':
                if statistics_3_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    from taxontabletools.create_krona_chart import create_krona_chart
                    create_krona_chart(statistics_3_taxon_table_path, path_to_outdirs)

            if event == 'run_CCA_analysis':

                if statistics_4_taxon_table_path == '':
                    sg.PopupError("Please provide a file")
                    print("Error: Please provide a file")

                else:
                    from taxontabletools.get_available_meta_data import get_available_meta_data
                    meta_data_to_test = get_available_meta_data(statistics_4_taxon_table_path, path_to_outdirs)

                    if meta_data_to_test == False:
                            sg.PopupError("Metadata file not found!")
                            print("Metadata file not found!")
                            print("\n" + "########")
                    else:
                            win2_active = True
                            window.Hide()
                            meta_data_to_test_list = list(slices([sg.Radio(name, "metadata", key=name, default=True) for name in sorted(meta_data_to_test)], 5))

                            layout2 = [[sg.Text("CCA analysis", size=(20,1))],
                                    [sg.Frame(layout = meta_data_to_test_list, title = 'Check taxa to exlude')],
                                    [sg.Button('Calculate')],
                                    [sg.Button('Back')]] # test = True

                            win2 = sg.Window('CCA analysis', layout2, keep_on_top=False)

                            while True:
                                event2, values2 = win2.Read()
                                if event2 == 'Calculate':
                                        for input_value, test in values2.items():
                                                if test == True:
                                                        meta_data_to_test = input_value
                                                        from taxontabletools.CCA_analysis import CCA_analysis
                                                        CCA_analysis(statistics_4_taxon_table_path, meta_data_to_test, cca_w, cca_h, cca_s, path_to_outdirs)
                                                        win2.Close()
                                                        win2_active = False
                                                        window.UnHide()
                                                        break

                                if event2 is None or event2 == 'Back':
                                        win2.Close()
                                        win2_active = False
                                        window.UnHide()
                                        break

            if event == 'run_PCoA_analysis':

                if statistics_4_taxon_table_path == '':
                    sg.PopupError("Please provide a file")
                    print("Error: Please provide a file")

                else:
                    from taxontabletools.get_available_meta_data import get_available_meta_data
                    meta_data_to_test = get_available_meta_data(statistics_4_taxon_table_path, path_to_outdirs)

                    if meta_data_to_test == False:
                            sg.PopupError("Metadata file not found!")
                            print("Metadata file not found!")
                            print("\n" + "########")
                    else:
                            win2_active = True
                            window.Hide()
                            meta_data_to_test_list = list(slices([sg.Radio(name, "metadata", key=name, default=True) for name in sorted(meta_data_to_test)], 5))

                            layout2 = [[sg.Text("PCoA analysis", size=(20,1))],
                                    [sg.Frame(layout = meta_data_to_test_list, title = 'Check taxa to exlude')],
                                    [sg.Button('Calculate')],
                                    [sg.Button('Back')]] # test = True

                            win2 = sg.Window('PCoA analysis', layout2, keep_on_top=False)

                            while True:
                                event2, values2 = win2.Read()
                                if event2 == 'Calculate':
                                        for input_value, test in values2.items():
                                                if test == True:
                                                        meta_data_to_test = input_value
                                                        from taxontabletools.PCoA_analysis import PCoA_analysis
                                                        PCoA_analysis(statistics_4_taxon_table_path, meta_data_to_test, pcoa_w, pcoa_h, pcoa_s, path_to_outdirs)
                                                        win2.Close()
                                                        win2_active = False
                                                        window.UnHide()
                                                        break

                                if event2 is None or event2 == 'Back':
                                        win2.Close()
                                        win2_active = False
                                        window.UnHide()
                                        break

            if event == 'run_replicate_consistency_filter':
                if taxon_table_processing_1_table == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    replicate_suffix_list = replicate_suffix_entry.replace(" ", "").split(",")
                    from taxontabletools.replicate_consistency_filter import replicate_consistency_filter
                    replicate_consistency_filter(taxon_table_processing_1_table, replicate_suffix_list, path_to_outdirs)

            if event == 'run_create_perlodes_input_file':
                if perlodes_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                    print("Error: Please provide a file")
                else:
                    from taxontabletools.convert_to_perlodes import convert_to_perlodes
                    convert_to_perlodes(perlodes_taxon_table_path, operational_taxon_list_path, path_to_outdirs)

            if event == 'run_alpha_diversity':

                    if statistics_4_taxon_table_path == '':
                        sg.PopupError("Please provide a file")
                        print("Error: Please provide a file")

                    else:
                        from taxontabletools.get_available_meta_data import get_available_meta_data
                        meta_data_to_test = get_available_meta_data(statistics_4_taxon_table_path, path_to_outdirs)

                        if meta_data_to_test == False:
                                sg.PopupError("Metadata file not found!")
                                print("Metadata file not found!")
                                print("\n" + "########")
                        else:
                                win2_active = True
                                window.Hide()
                                meta_data_to_test_list = list(slices([sg.Radio(name, "metadata", key=name, default=True) for name in sorted(meta_data_to_test)], 5))

                                layout2 = [[sg.Text("Alpha diversity", size=(20,1))],
                                [sg.Frame(layout = meta_data_to_test_list, title = 'Check meta data to test')],
                                [sg.Button('Calculate')],
                                [sg.Button('Back')]]

                                win2 = sg.Window('Alpha diversity estimates', layout2, keep_on_top=False)

                                while True:
                                    event2, values2 = win2.Read()
                                    if event2 == 'Calculate':
                                            for meta_data_to_test, test in values2.items():
                                                if test == True:
                                                    if alpha_box_plot == True:
                                                            from taxontabletools.alpha_diversity import alpha_diversity_boxplot
                                                            alpha_diversity_boxplot(statistics_4_taxon_table_path, meta_data_to_test, alpha_w, alpha_h, alpha_s, alpha_font, path_to_outdirs)
                                                            win2.Close()
                                                            win2_active = False
                                                            window.UnHide()
                                                            break
                                                    if alpha_scatter_plot == True:
                                                            from taxontabletools.alpha_diversity import alpha_diversity_scatter_plot
                                                            alpha_diversity_scatter_plot(statistics_4_taxon_table_path, meta_data_to_test, alpha_w, alpha_h, alpha_s, alpha_font, path_to_outdirs)
                                                            win2.Close()
                                                            win2_active = False
                                                            window.UnHide()
                                                            break
                                    if event2 is None or event2 == 'Back':
                                            win2.Close()
                                            win2_active = False
                                            window.UnHide()
                                            break

            if event == 'run_beta_diversity':

                    if statistics_4_taxon_table_path == '':
                        sg.PopupError("Please provide a file")
                        print("Error: Please provide a file")

                    else:
                        from taxontabletools.get_available_meta_data import get_available_meta_data
                        meta_data_to_test = get_available_meta_data(statistics_4_taxon_table_path, path_to_outdirs)

                        if meta_data_to_test == False:
                                sg.PopupError("Metadata file not found!")
                                print("Metadata file not found!")
                                print("\n" + "########")
                        else:
                                win2_active = True
                                window.Hide()
                                meta_data_to_test_list = list(slices([sg.Radio(name, "metadata", key=name, default=True) for name in sorted(meta_data_to_test)], 5))

                                layout2 = [[sg.Text("Beta diversity", size=(20,1))],
                                [sg.Frame(layout = meta_data_to_test_list, title = 'Check meta data to test')],
                                [sg.Button('Calculate')],
                                [sg.Button('Back')]]

                                win2 = sg.Window('Beta diversity estimates', layout2, keep_on_top=False)

                                while True:
                                    event2, values2 = win2.Read()
                                    if event2 == 'Calculate':
                                            for meta_data_to_test, test in values2.items():
                                                if test == True:
                                                    from taxontabletools.beta_diversity import beta_diversity
                                                    beta_diversity(statistics_4_taxon_table_path, beta_w, beta_h, beta_font, beta_cmap, meta_data_to_test, path_to_outdirs)
                                                    win2.Close()
                                                    win2_active = False
                                                    window.UnHide()
                                                    break
                                    if event2 is None or event2 == 'Back':
                                            win2.Close()
                                            win2_active = False
                                            window.UnHide()
                                            break

            if event == 'run_log_network':
                from taxontabletools.create_log import ttt_log_network
                ttt_log_network(path_to_outdirs)

            if event == 'open_fgbewertung':
                print("Open: https://www.gewaesser-bewertung-berechnung.de")
                webbrowser.open('https://www.gewaesser-bewertung-berechnung.de')

            if event == 'open_operational_taxonlist_download':
                print("Open: https://www.gewaesser-bewertung-berechnung.de/index.php/perlodes-online.html")
                webbrowser.open('https://www.gewaesser-bewertung-berechnung.de/index.php/perlodes-online.html')

            if event == 'open_github':
                print("Open: https://github.com/TillMacher/TaxonTableTools")
                webbrowser.open('https://github.com/TillMacher/TaxonTableTools')

            if event == 'open_twitter':
                print("Open: https://twitter.com/THM_93")
                webbrowser.open('https://twitter.com/THM_93')

            if event == 'run_metadata_table_help':
                sg.Popup(metadata_table_help_text(), title = "Metadata table")

            if event == 'run_alpha_diversity_help_text':
                sg.Popup(alpha_diversity_help_text(), title = "Alpha diversity")

            if event == 'run_beta_diversity_help_text':
                sg.Popup(beta_diversity_help_text(), title = "Beta diversity")

            if event == 'run_cca_analysis_help_text':
                sg.Popup(cca_analysis_help_text(), title = "CCA analysis")

            if event == 'run_pcoa_analysis_help_text':
                sg.Popup(pcoa_analysis_help_text(), title = "PCoA analysis")

            # close all windows to prevent crazy overlaps
            if plt.get_fignums():
                plt.close()

            print("\n" + "########")

        ###############################################################################################
        # exception collection area

        # end the script when KeyboardInterrupt is raised
        except KeyboardInterrupt:
            exit()

        # print error message if a module is not installed
        except ModuleNotFoundError:
            if win2_active == True:
                win2.Close()
                win2_active = False
                window.UnHide()

            sg.PopupError("Error: Please install all required modules.", keep_on_top=True)
            print("\n" + "########")

        # Raised when an error is detected that doesn’t fall in any of the other categories. The associated value is a string indicating what precisely went wrong.
        # CHANGE TO CUSTOM ERROR and not RuntimeError!!!
        except RuntimeError:
            print("\n" + "########")
            try:
                window.UnHide()
            except:
                pass
            try:
                win2.Close()
            except:
                pass


        # if there are unexpected errors print a message and continue the script!
        except:
            print("\n" + "########")
            # UnHide the main window if neccessary
            try:
                window.UnHide()
            except:
                pass
            try:
                win2.Close()
            except:
                pass

            layout = [
                        [sg.Image(crash), sg.Text(" You've been awarded with the gold medal in program crashing!")],
                        [sg.Text("", size=(1,1))],
                        [sg.Text("Unexpected error: " + str(sys.exc_info()[0]))],
                        [sg.Text("", size=(1,1))],
                        [sg.Text("An unexpected error occured:")],
                        [sg.Text("Please refer to the manual.")],
                        [sg.Text("", size=(1,1))],
                        [sg.Button('Return'), sg.Button('Exit', button_color=('black', 'red'))]
                     ]

            window_error = sg.Window('Error', layout, keep_on_top=True)

            while (True):

                event, values = window_error.Read()

                if event in ('Return', None):
                    break

                if event in ('Exit', None):
                    raise

            window_error.Close()

            print("\n" + "########")

        ###############################################################################################

    window.Close()

## run only if called as toplevel script
if __name__ == "__main__":
    main()
