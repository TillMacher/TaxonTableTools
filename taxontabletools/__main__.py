import glob, sys, os, pkgutil, ast, subprocess
import PySimpleGUI as sg
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
from pathlib import Path
import webbrowser
from lastversion import lastversion
import feedparser

##########################################################################################################################
# update version here (will be displayed on the main layout)
taxon_tools_version = str('1.3.1')

##########################################################################################################################
# general functions

# slice function for lists to split up lists
def slices(list, slice):
    for i in range(0, len(list), slice):
        yield list[i : i + slice]

def open_taxon_table(taxon_table):
    if sys.platform == "win32":
        os.startfile(taxon_table)
    else:
        opener = "open" if sys.platform == 'darwin' else 'xdg-open'
        subprocess.call([opener, taxon_table])

##########################################################################################################################
##########################################################################################################################

def change_log_text():
    change_log_text = """
    TTT change log
    ____________________________________________________________________________

    v 1.3.1
    New functions:
    - Metadata entries for specific samples can now be left blank to exclude the
      sample from the analyses. This avoids the need to create new metadata tables
      for each subset of a TaXon table
    - Taxon lists have been overhauled and now include presence/absence indications
      for samples or metadata  groups
    Bug fixes:
    - Fixed incorrect y-axis title in the per taxon plot
    - Enhanced y-axis title of shared OTUs per bin plot
    - Fixed Qiime2 conversion
    - Fixed OTU and read correlation layout bugs

    """
    return change_log_text

def metadata_table_help_text():
    metadata_table_help_text = """
    A metadata will be automatically created as Excel sheet.
    ____________________________________________________________________________

    Please rename the columns according to your metadata in the Excel sheet!
    This cannot be performed via TTT! The default column names are non-descriptive
    and function as placeholders!

    The columns will be read from the table and will show up in TTT to select
    the metadata to test on.

    Each value should be adjusted to the users need before loading the metadata table.

    """
    return metadata_table_help_text

def alpha_diversity_help_text():
    alpha_diversity_help_text = """

    Alpha diversity is calculated with the python package scikit-bio (http://scikit-bio.org/)
    and is based on the number of taxa (user choice) per sample, which are displayed as scatter or bar plot.

    """
    return alpha_diversity_help_text

def beta_diversity_help_text():
    beta_diversity_help_text = """

    Beta diversity is calculated with the python package scikit-bio (http://scikit-bio.org/)
    and is either calculated as jaccard or braycurtis distances, which are illustrated in a distance matrix.
    An R-square test for reliability and validity is conducted.

    Jaccard distances:
    Based on boolean (presence/absence data) >> Recommended

    Braycurtis distances:
    Based on read counts >> Not recommended. Reads do not necessarily reflect specimen numbers or read counts.
    Consider this when interpreting the results.

    The TaXon table data will be automatically converted to incidence (boolean) values.

    """
    return beta_diversity_help_text

def cluster_dendrogram_help_text():
    cluster_dendrogram_help_text = """

    The cluster dendrogram is calculated with the python package scikit-bio (http://scikit-bio.org/)
    and is based on jaccard or braycurtis distances.
    The group threshold can be used to highlight groupings in the dendrogram.

    Jaccard distances:
    Based on boolean (presence/absence data) >> Recommended

    Braycurtis distances:
    Based on read counts >> Not recommended. Reads do not necessarily reflect specimen numbers or read counts.
    Consider this when interpreting the results.

    The TaXon table data will be automatically converted to incidence (boolean) values.

    """
    return cluster_dendrogram_help_text

def cca_analysis_help_text():
    cca_analysis_help_text = """

    Canonical-correlation analysis (CCA) analysis plots are calculated with the python package scikit-bio (http://scikit-bio.org/).
    The two axes to plot can be chosen from all available axes.
    The metadata must be a vector of random variables (X = X1...,Xn) and must differ between samples!
    All eigenvalues are saved as an Excel sheet.

    The TaXon table data will be automatically converted to incidence (boolean) values.

    """
    return cca_analysis_help_text

def pcoa_analysis_help_text():
    pcoa_analysis_help_text = """

    Principle coordinate analysis (PCoA) plots are calculated with the python package scikit-bio (http://scikit-bio.org/)
    and are based on jaccard or braycurtis distances.
    The metadata must differ between samples!
    An R-square test for reliability and validity is conducted.
    The two axes to plot can be chosen from all available axes.
    All eigenvalues are saved as an Excel sheet.

    Jaccard distances:
    Based on boolean (presence/absence data) >> Recommended

    Braycurtis distances:
    Based on read counts >> Not recommended. Reads do not necessarily reflect specimen numbers or read counts.
    Consider this when interpreting the results.

    The TaXon table data will be automatically converted to incidence (boolean) values.

    """
    return pcoa_analysis_help_text

def nmds_analysis_help_text():
    nmds_analysis_help_text = """

    Non-metric multidimensional scaling analysis (NMDS) plots are calculated with the python package scikit-bio (http://scikit-bio.org/)
    and are based on jaccard-distances. The TaXon table data will be automatically converted to incidence (boolean) values.
    The metadata must differ between samples!

    A stress value will be calculated to measure the reliability of the results:
        -> Stress values above 0.3 are considered arbitrary.
        -> Stress values between 0.2 and 0.1 are considered suspect.
        -> Stress values below 0.1 are considered to indicate reliable patterns.
        -> Stress values below 0.05 are considered a good fit.

    If you experience high stress values, increase the number of initiations and maximum number of iterations.

    From the scikit manual:
    Number of initiations:
    'Number of times the SMACOF algorithm will be run with different initializations. The final results will be the best output of the runs,
    determined by the run with the smallest final stress.'
    Maximum number of iterations:
    'Maximum number of iterations of the SMACOF algorithm for a single run.'

    Jaccard distances:
    Based on boolean (presence/absence data) >> Recommended
    Braycurtis distances:
    Based on read counts >> Not recommended. Reads do not necessarily reflect specimen numbers or read counts.
    Consider this when interpreting the results.

    """
    return nmds_analysis_help_text

def gbif_check_help_text():
    gbif_check_help_text = """

    The assigned taxonomy of each OTU in the TaXon table will be checked against the GBIF (Global Biodiversity Information Facility) database,
    using the API.

    This tool will recognize spelling mistakes and search for synonyms for the respective taxon. Changes will be automatically be applied
    in the TaXon table. A respective log file will be saved to track the changes in the taxonomy.

    Disclaimer: The GBIF database is a highly complete and curated database, but mistakes can slip in during the taxonomy matching.
    Thus, please check the changes in detail and adjust wrong assignments manually.

    """
    return gbif_check_help_text

def pa_conversion_help_text():
    pa_conversion_help_text = """
    The TaXon table will be converted to incidence data.
    This format is required for most diversity and ordiniation analyses.
    ____________________________________________________________________________

    The use of read abundances as a proxy for specimen counts or biomass estimates
    has been subject of discussion with the development of DNA metabarcoding.
    Due to PCR stochasticity, varying primer binding efficiency and sequencing bias,
    there is often only a weak correlation between read abundances and specimen
    counts or biomass, although for several cases with improved primer settings
    exceptions exist. Thus, often it is recommended to convert the read abundance
    data to incidence data for biodiversity analyses. However, this conversion
    comes with a downside, since incidence data limits the pool of appropriate
    diversity estimate analyses.

    """
    return pa_conversion_help_text

def table_per_sample_help_text():
    table_per_sample_help_text = """
    Each sample from the selected TaXon table will be assigned to its own respective table
    ____________________________________________________________________________

    This tool is particularly useful to compare single samples in the Venn diagram tool.
    OTUs that are not present in the respective sample are discarded.

    """
    return table_per_sample_help_text

def read_proportions_help_text():
    read_proportions_help_text = """
    Read proportion plots, although informative, can quickly become overwhelming, difficult to read and hard to interprete.
    ____________________________________________________________________________

    The easiest workaround is to reduce the amount of data that is presented. Follow these
    few steps and the read proportions plots will quickly become much more eye-pleasing!
    1. Choosing a higher taxonomic taxonomic level
        -> On lower taxonomic levels (e.g. species) too many taxa are presented at once!
    2. Filter the data
        -> Removal of not required samples
        -> Removal of the groups out of the study's scope
        -> Keep only species level hits (see Tutorial how to do that)

    """
    return read_proportions_help_text

def layout_help_text():
    layout_help_text = """

    TTT offers a wide range of customizable functions to individualize plots.
    ____________________________________________________________________________

    Layout:
    A layout defines the general appearance of a plot, such as the grid,
    the background or axes ticks. Various widely distributed layouts are available,
    for example 'seaborn' and 'ggplot2'.

    Themes:
    A theme defines the two primary colors of a plot. Several pre-defined
    themes are available. Additionally, a custom theme consisting of two colors
    can be manually selected.

    Color sequences (qualitative):
    A color sequence defines a range of discrete colors, that is for example used in bar charts
    or scatter plots.

    Color scales (sequential):
    A color scale defines a range of sequential colors, that is for example used in heatmaps.

    """
    return layout_help_text

##########################################################################################################################
##########################################################################################################################

def main():

    sg.ChangeLookAndFeel("Reddit")

    ## check for the latest version of TTT
    try:
        latest_version = lastversion.has_update(repo="taxontabletools", at='pip', current_version=taxon_tools_version)
        if latest_version:
            sg.Popup('Newer TTT version is available: ' + str(latest_version) + "\nPlease update TTT via pip!", title="Update available!")
    except:
        pass

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
    fliessgewaesserbewertungde = Path(str(path_to_ttt) + "/_source/fliessgewaesserbewertung_de.png")
    github = Path(str(path_to_ttt) + "/_source/github.png")
    ttt_main = Path(str(path_to_ttt) + "/_source/main.png")
    ttt_logo = Path(str(path_to_ttt) + "/_source/taxontabletools.png")
    twitter = Path(str(path_to_ttt) + "/_source/twitter.png")
    help_button = Path(str(path_to_ttt) + "/_source/help_button.png")

    # fresh start: there is an empty user_data file
    # ask for user Input
    # stay open until a path was defined
    # then write it the user_data file to reload
    while projects_main_path == "":
        sg.Popup(change_log_text(), title = "TTT change log")
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

    ## load the user preferences
    try:
        user_preferences = Path(str(path_to_ttt) + "/user_preferences.txt")
        f = open(user_preferences, "r")
        default_template = f.readline().strip()
        default_color1 = f.readline().strip()
        default_color2 = f.readline().strip()
        default_color_sequence = f.readline().strip()
        default_colorscale = f.readline().strip()
        default_font_size = int(f.readline().strip())
        default_transparent_theme = f.readline().strip()
        default_clustering_unit = f.readline().strip()
        if default_transparent_theme == "True":
            default_transparent_theme = True
        else:
            default_transparent_theme = False
        f.close()
    except:
        default_template = "simple_white"
        default_theme = "Blue"
        default_color1 = "Blue"
        default_color2 = "Black"
        default_color_sequence = "Plotly"
        default_colorscale = "viridis"
        default_font_size = 16
        default_transparent_theme = False
        default_clustering_unit = "OTUs"

    taXon_table_dir_path = Path(str(path_to_outdirs) + "/TaXon_tables/")

    directories_to_create = ["Venn_diagrams","TaXon_tables", "TaXon_tables_per_sample", "Taxon_lists", "Rarefaction_curves",
    "Taxonomic_richness_plots", "CCA_plots", "Taxonomic_resolution_plots", "Meta_data_table", "Site_occupancy_plots", "Read_proportions_plots", "Basic_stats",
    "Krona_charts", "Perlodes", "Alpha_diversity", "Beta_diversity", "PCoA_plots", "Replicate_analysis", "GBIF", "Occurrence_analysis", "Per_taxon_statistics",
    "ParCat_plots", "NMDS_plots"]

    for directory in directories_to_create:
        dirName = Path(str(path_to_outdirs) + "/" + directory + "/")
        if not os.path.exists(dirName):
            os.mkdir(Path(dirName))

    available_templates_list = ['seaborn', 'ggplot2', 'simple_white', 'plotly', 'plotly_dark', 'presentation', 'plotly_white']
    available_color1_list = ["Dark", "Light", "Green", "Blue"]

    plotly_colors = ["aliceblue", "antiquewhite", "aqua", "aquamarine", "azure",
    "beige", "bisque", "black", "blanchedalmond", "blue",
    "blueviolet", "brown", "burlywood", "cadetblue",
    "chartreuse", "chocolate", "coral", "cornflowerblue",
    "cornsilk", "crimson", "cyan", "darkblue", "darkcyan",
    "darkgoldenrod", "darkgray", "darkgrey", "darkgreen",
    "darkkhaki", "darkmagenta", "darkolivegreen", "darkorange",
    "darkorchid", "darkred", "darksalmon", "darkseagreen",
    "darkslateblue", "darkslategray", "darkslategrey",
    "darkturquoise", "darkviolet", "deeppink", "deepskyblue",
    "dimgray", "dimgrey", "dodgerblue", "firebrick",
    "floralwhite", "forestgreen", "fuchsia", "gainsboro",
    "ghostwhite", "gold", "goldenrod", "gray", "grey", "green",
    "greenyellow", "honeydew", "hotpink", "indianred", "indigo",
    "ivory", "khaki", "lavender", "lavenderblush", "lawngreen",
    "lemonchiffon", "lightblue", "lightcoral", "lightcyan",
    "lightgoldenrodyellow", "lightgray", "lightgrey",
    "lightgreen", "lightpink", "lightsalmon", "lightseagreen",
    "lightskyblue", "lightslategray", "lightslategrey",
    "lightsteelblue", "lightyellow", "lime", "limegreen",
    "linen", "magenta", "maroon", "mediumaquamarine",
    "mediumblue", "mediumorchid", "mediumpurple",
    "mediumseagreen", "mediumslateblue", "mediumspringgreen",
    "mediumturquoise", "mediumvioletred", "midnightblue",
    "mintcream", "mistyrose", "moccasin", "navajowhite", "navy",
    "oldlace", "olive", "olivedrab", "orange", "orangered",
    "orchid", "palegoldenrod", "palegreen", "paleturquoise",
    "palevioletred", "papayawhip", "peachpuff", "peru", "pink",
    "plum", "powderblue", "purple", "red", "rosybrown",
    "royalblue", "saddlebrown", "salmon", "sandybrown",
    "seagreen", "seashell", "sienna", "silver", "skyblue",
    "slateblue", "slategray", "slategrey", "snow", "springgreen",
    "steelblue", "tan", "teal", "thistle", "tomato", "turquoise",
    "violet", "wheat", "white", "whitesmoke", "yellow",
    "yellowgreen"]

    available_colorsequences = {"Plotly":px.colors.qualitative.Plotly, "G10":px.colors.qualitative.G10,
    "T10":px.colors.qualitative.T10, "Alphabet":px.colors.qualitative.Alphabet, "Dark24":px.colors.qualitative.Dark24, "Dark24_r":px.colors.qualitative.Dark24_r,
    "Light24":px.colors.qualitative.Light24, "Set1":px.colors.qualitative.Set1, "Pastel1":px.colors.qualitative.Pastel,
    "Dark2":px.colors.qualitative.Dark2, "Set2":px.colors.qualitative.Set2, "Pastel2":px.colors.qualitative.Pastel2,
    "Set3":px.colors.qualitative.Set3, "Antique":px.colors.qualitative.Antique,"Bold":px.colors.qualitative.Bold,
    "Pastel":px.colors.qualitative.Pastel, "Prism":px.colors.qualitative.Prism, "Safe":px.colors.qualitative.Safe,
    "Vivid":px.colors.qualitative.Vivid}

    available_colorscales = px.colors.named_colorscales()

    available_taxonomic_levels_list = ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'OTUs']

    ##########################################################################################################################
    ##########################################################################################################################
    ##########################################################################################################################

    ## define colum layouts for the layout options
    ## column1
    column1 = [[sg.Text('Layout:')], [sg.Combo(available_templates_list, default_value=default_template, key="template", size=(10,1))]]
    ## column2
    column2 = [[sg.Text('Theme:')], [sg.Combo(plotly_colors, default_value=default_color1, key="color1", size=(12,1))], [sg.Combo(plotly_colors, default_value=default_color2, key="color2", size=(12,1))]]
    ## column3
    column3 = [[sg.Text("Color sequence:")], [sg.Combo(list(available_colorsequences.keys()), default_value=default_color_sequence, key="colorsequence", size=(12,1))],
                [sg.Button("Preview", key = 'run_show_colorsequence', button_color=('black', 'white'))],
                [sg.Text('')],
                [sg.Text("Color scale:")], [sg.Combo(available_colorscales, default_value=default_colorscale, key="colorscale", size=(12,1))],
                [sg.Button("Preview", key = 'run_show_colorscale', button_color=('black', 'white'))]]
    ##column4
    OTU_IDs = ["ASVs", "ESVs", "OTUs", "zOTUs"]
    column4 = [[sg.Text('Options:')], [sg.Text("Transparent theme:", size=(15,1)), sg.CB("", key="transparent_theme", default=default_transparent_theme)],
                [sg.Text("Font size:", size=(15,1)), sg.Combo([i for i in range(6, 22, 2)], default_value=default_font_size, key="font_size", size=(12,1))],
                [sg.Text('Clustering unit:', size=(15,1)), sg.Combo(OTU_IDs, default_value=default_clustering_unit, key="clustering_unit", size=(12,1))]]
    ##column5
    column5 = [[sg.Text('')],
                [sg.Button(" ? ", key = 'run_layout_help', button_color=('blue', 'white'))],
                [sg.Button("Save preferences", key = 'run_save_preferences', button_color=('black', 'white'))]]

    getting_started_layout = [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('TaXon table converter', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('Enter path to taxonomy table and read table to create a new table in TaXon format', size=(85,2))],
    					[sg.Text('Taxonomy table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'taxonomy_table_path'), sg.Button("Check table", key = 'run_taxonomy_table_check', button_color=('black', 'white'))],
                        [sg.Text('Sorting method:', size=(20, 1)),  sg.Radio('First hit', "sorting_method", key='sorting_method_fh'), sg.Radio('JAMP', "sorting_method", default=True, key='sorting_method_jamp'), sg.Radio('BOLDigger', "sorting_method", key='sorting_method_boldigger')],
    					[sg.Text('Read table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'Read_table_path'), sg.Button("Check table", key = 'run_read_table_check', button_color=('black', 'white'))],
                        [sg.Text('Format:', size=(20, 1)), sg.Radio('TTT', "read_table_format", key='read_table_format_TTT', default=True),  sg.Radio('Qiime2', "read_table_format", key='read_table_format_qiime')],
                        [sg.Text('Name of TaXon table:', size=(20, 1)), sg.Input(project_folder + "_taxon_table", key = 'TaXon_table_name'), sg.Text('.xlsx'), sg.Button("Convert", key = 'run_taxon_converter')],
    					[sg.Text('')],
                        [sg.Text('Check TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(initial_folder = taXon_table_dir_path, key = 'taxon_table_0_path'), sg.Button("Check table", key = 'run_taXon_table_check', button_color=('black', 'white')), sg.Button("Open", key = 'open0')],
                        [sg.Text('_'*120)],
    					[sg.Text('')],
    					[sg.Text('Graphic layout options', size=(50,2), font=('Arial', 11, "bold"))],
                        [sg.Column(column1), sg.Column(column2), sg.Column(column3), sg.Column(column4), sg.Column(column5)],
    					[sg.Text('')],
    					]

    replicate_analyses_layout = [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Replicate analyses', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(initial_folder = taXon_table_dir_path, key = 'taxon_table_1_path'), sg.Button("Open", key="open1")],
    					[sg.Text('_'*120)],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Replicates',size=(40,2), font=('Arial', 11, "bold"))],
                        [sg.Text("Replicate suffixes:", size=(30, 1)), sg.Input("a, b", size=(15,1), key="replicate_suffix_entry"),
                        sg.Text("", size=(1,1)), sg.Button("Preview", key="run_preview_suffix"), sg.Text("", size=(1,1)),
                        sg.Frame(layout=[[sg.Text("Comma separated list of suffixes.")]], title="Information")],
                        [sg.Text('',size=(1,1))],
    					[sg.Text('Replicate analysis',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Calculate shared OTUs", size=(30,1)), sg.Button("Run", key = 'run_replicate_analysis'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("1000", size=(4,1), key="width_repicates"),
                        sg.Input("600", size=(4,1), key="height_repicates"), sg.CB("Custom colors:", default=False, key="replicate_analysis_custom_colors")]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text("Calculate OTU and read correlation", size=(30,1)), sg.Button("Run", key = 'run_replicate_correlation'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("1200", size=(4,1), key="width_repcorr"),
                        sg.Input("600", size=(4,1), key="height_repcorr"), sg.Text("", size=(1,1)), sg.Text("Start axes at zero"), sg.CB("x", key="x_zero_repcorr"), sg.CB("y", key="y_zero_repcorr")]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Replicate processing',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Merge replicates", size=(30, 1)), sg.Button("Run", key='run_replicate_consistency_filter'), sg.Text("", size=(1,1)),
                        sg.Frame(layout=[[sg.Text("Merge replicates and optionally keep only\nOTUs that are present in all replicates."), sg.CB("Replicate consistency", key = "consistency_option", default=True)]], title="Information")],
    					[sg.Text('',size=(1,1))],
    					]

    taxon_table_filtering_layout = [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('TaXon table filtering', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(initial_folder = taXon_table_dir_path, key = 'taxon_table_2_path'), sg.Button("Open", key="open2")],
    					[sg.Text('_'*120)],
    					[sg.Text('',size=(1,1))],

    					[sg.Text('Taxon-based filter',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Taxon-based filter:", size=(20, 1)), sg.Button("Run", key = 'run_taxon_filtering'), sg.Text("", size=(2,1)), sg.Frame(layout=[[sg.Radio('all', "cm_taxon_based_filter", key='cm_taxon_based_filter_all', default=True),
                        sg.Radio('none', "cm_taxon_based_filter", key='cm_taxon_based_filter_none')]], title="Check marks"), sg.Frame(layout=[[sg.Radio('Keep', "taxon_filter_method", key='taxon_filter_method_keep'),
                        sg.Radio('Exclude', "taxon_filter_method", key='taxon_filter_method_exclude', default=True)]], title="Method"), sg.Text("Mask:"), sg.Input(size=(5,1), key="taxon_filter_method_criterion"),
                        sg.Text('Tax. level:'), sg.Combo(available_taxonomic_levels_list, default_value='Phylum', key="filter_taxonomic_level")],

    					[sg.Text('',size=(1,1))],

    					[sg.Text('Sample-based filter', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Sample-based filter:", size=(20, 1)), sg.Button("Run", key = 'run_sample_filtering'), sg.Text("", size=(2,1)), sg.Frame(layout=[[sg.Radio('all', "cm_sample_based_filter", key='cm_sample_based_filter_all', default=True),
                        sg.Radio('none', "cm_sample_based_filter", key='cm_sample_based_filter_none')]], title="Check marks"), sg.Frame(layout=[[sg.Radio('Keep', "sample_filter_method", key='sample_filter_method_keep'),
                        sg.Radio('Exclude', "sample_filter_method", key='sample_filter_method_exclude', default=True)]], title="Method"), sg.Text("Mask:"), sg.Input(size=(5,1), key="sample_filter_method_criterion")],

    					[sg.Text('',size=(1,1))],

    					[sg.Text('Read-based filter', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('Read-based filter:', size=(20, 1)), sg.Button("Run", key='run_read_filter'), sg.Text("", size=(1,1)),
                        sg.Frame(layout=[[sg.Text('Treshold:'), sg.Radio('Relative', "read_filtering_radio", key="relative_filtering", default=True), sg.Input('0.01', size=(5,1), key="relative_filter_threshold"), sg.Text("%"),
                        sg.Radio('Absolute', "read_filtering_radio", key="absolute_filtering"), sg.Input('50', size=(5,1), key="absolute_filter_threshold")]], title="Settings")],

    					[sg.Text('',size=(1,1))],

    					[sg.Text('Negative control subtraction', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Subtract NC reads :", size=(20, 1)), sg.Button("Run", key = 'run_nc_subtraction'), sg.Text("", size=(2,1)), sg.Frame(layout=[[sg.Radio('all', "cm_nc_subtraction", key='cm_nc_subtraction_all', default=True),
                        sg.Radio('none', "cm_nc_subtraction", key='cm_nc_subtraction_none')]], title="Check marks"), sg.Text("Search for negative controls:"), sg.Input(size=(7,1), key="nc_subtraction_criterion")],

    					]

    data_conversion_layout = [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Data conversion', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(initial_folder = taXon_table_dir_path, key = 'taxon_table_3_path'), sg.Button("Open", key="open3")],
    					[sg.Text('_'*120)],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('TaXon table per sample', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Create TaXon tables per-sample", size=(25,1)), sg.Button("Run", key = 'run_taxon_table_per_sample'), sg.Text("", size=(1,1)),
                        sg.Button(" ? ", button_color=('blue', 'white'), key="run_table_per_sample_help" )],
    					[sg.Text('',size=(1,1))],
                		[sg.Text('Presence absence table', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Convert to p/a data", size=(25,1)), sg.Button("Run", key = 'run_convert_to_pa'), sg.Text("", size=(1,1)),
                        sg.Button(" ? ", button_color=('blue', 'white'), key="run_pa_conversion_help" )],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Metadata table', size=(50,2), font=('Arial', 11, "bold")), sg.Text("", size=(2,1))],
    					[sg.Text("Create metadata table", size=(25,1)), sg.Button("Run", key = 'run_create_meta_data_table'), sg.Text("", size=(1,1)),
                        sg.Button(" ? ", button_color=('blue', 'white'), key="run_metadata_table_help" ), sg.Button("Modify", button_color=('blue', 'white'), key = 'run_modify_metadata_table')],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('GBIF taxonomy check', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Check TaXon table:", size=(25, 1)), sg.Button("Run", key = 'run_gbif_check_taxonomy'), sg.Text("", size=(1,1)),
                        sg.Button(" ? ", button_color=('blue', 'white'), key="run_gbif_check_help" )],
    					[sg.Text('',size=(1,1))],
    					]

    basic_statistics_layout = [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Basic statistics', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(initial_folder = taXon_table_dir_path, key = 'taxon_table_4_path'), sg.Button("Open", key="open4")],
    					[sg.Text('_'*120)],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Basic statistics',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Calculate basic statistics", size=(25,1)), sg.Button("Run", key = 'run_basic_stats'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("800", size=(4,1), key="width_bstats"),
                        sg.Input("1100", size=(4,1), key="height_bstats"), sg.Text("", size=(1,1))]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Taxonomic resolution',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Taxonomic resolution plot", size=(25,1)), sg.Button("Run", key = 'run_taxonomic_resolution'), sg.Text("", size=(1,1)),
                        sg.Frame(layout=[
                        [sg.Text("Plot size (w,h):"), sg.Input("600", size=(4,1), key="width_tax_res"), sg.Input("600", size=(4,1), key="height_tax_res"), sg.Text("Plot type:"), sg.Radio('a', "tres_plot_type", key='tres_type_a', default=True), sg.Radio('b', "tres_plot_type", key='tres_type_b')]
                        ], title="Settings", size=(1,2))],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Taxonomic richness',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Taxonomic richness plot", size=(25,1)), sg.Button("Run", key = 'run_taxonomic_richness'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("600", size=(4,1), key="width_tax_rich"),
                        sg.Input("600", size=(4,1), key="height_tax_rich")]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					]

    sample_comparison_layout =  [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Sample comparison', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(initial_folder = taXon_table_dir_path, key = 'taxon_table_5_path'), sg.Button("Open", key="open5")],
    					[sg.Text('_'*120)],
    					[sg.Text('',size=(1,1))],

    					[sg.Text('Venn diagrams',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Compare two TaXon tables:", size=(25,1)), sg.Button("Run", key = 'run_venn2_diagram'), sg.Text("", size=(2,1)), sg.Text("2nd TaXon table:"), sg.Input("", size=(15,1)), sg.FileBrowse(initial_folder = taXon_table_dir_path, key = 'venn_taxon_table_2_path')],
                        [sg.Text("Compare three TaXon tables:", size=(25,1)), sg.Button("Run", key = 'run_venn3_diagram'), sg.Text("", size=(2,1)), sg.Text("3rd TaXon table:"), sg.Input("", size=(15,1)), sg.FileBrowse(initial_folder = taXon_table_dir_path, key = 'venn_taxon_table_3_path')],
    					[sg.Text('',size=(1,1))],
                        [sg.Text("Output file: "), sg.Input("venn_diagram", size=(15,1), key="venn_diagram_name")],

    					[sg.Text('',size=(1,1))],
    					[sg.Text('Rarefaction curve',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Sample-based rarefaction", size=(25,1)), sg.Button("Run", key = 'run_rarefaction_curve'), sg.Text("", size=(1,1)),
                        sg.Frame(layout=[[sg.Text("Method:"), sg.Combo(["all-in-one", "per taxon"], default_value="all-in-one", key="rarefaction_curve_method"), sg.Text('Repetitions:'), sg.Input('1000', size=(6,1), key='sample_repetitions')],
                        [sg.Text('Tax. level to test:'), sg.Combo(available_taxonomic_levels_list, default_value='Species', key="rarefaction_curve_taxonomic_level_1"),
                        sg.Text('Tax. level to separate:'), sg.Combo(available_taxonomic_levels_list, default_value='Class', key="rarefaction_curve_taxonomic_level_2")]], title="Settings")],

                        [sg.Text("Read-based rarefaction", size=(25,1)), sg.Button("Run", key = 'run_rarefaction_curve_reads'), sg.Text("", size=(1,1)),
                        sg.Frame(layout=[[sg.Text('Repetitions:'), sg.Input('10', size=(6,1), key='read_repetitions'), sg.Text("Plot size (w,h): "), sg.Input("1000", size=(4,1), key="width_rarefaction_reads"),
                        sg.Input("1000", size=(4,1), key="height_rarefaction_reads"), sg.Text("", size=(1,1))]], title="Settings")],

    					[sg.Text('',size=(1,1))],
    					[sg.Text('Site occupancy', size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Calculate site occupancy", size=(25,1)), sg.Button("Run", key = 'run_site_occupancy'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h): "), sg.Input("1000", size=(4,1), key="width_site_occ"),
                        sg.Input("1000", size=(4,1), key="height_site_occ"), sg.Text("Tax. level:"), sg.Combo(available_taxonomic_levels_list, default_value="Species", key="site_occupancy_taxonomic_level")],
                        [sg.Text("Layout:"), sg.Combo(["barchart", "heatmap"], default_value="heatmap", key="site_occupancy_option"), sg.CB("Summarise metadata", default=True, key="site_occupancy_heatmap_add_metadata_sum")]], title="Settings")],

    					[sg.Text('',size=(1,1))],
    					]

    taxonomic_composition_layout =  [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Taxonomic composition', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(initial_folder = taXon_table_dir_path, key = 'taxon_table_6_path'), sg.Button("Open", key="open6")],
    					[sg.Text('_'*120)],

    					[sg.Text('',size=(1,1))],
    					[sg.Text('Per taxon statistics',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Reads, OTUs and species", size=(25,1)), sg.Button("Run", key = 'run_per_taxon_analysis'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("1400", size=(4,1), key="width_per_taxon_stats"),
                        sg.Input("700", size=(4,1), key="height_per_taxon_stats"), sg.Text('Tax. level:'), sg.Combo(['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'], default_value="Class", key="per_taxon_stats_taxonomic_level")]], title="Settings")],

    					[sg.Text('',size=(1,1))],
    					[sg.Text('Parallel category analysis',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Calculate ParCat plot", size=(25,1)), sg.Button("Run", key = 'run_parcat_analysis'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("1000", size=(4,1), key="width_parcat"),
                        sg.Input("1200", size=(4,1), key="height_parcat")],
                        [sg.Text('Tax. level:'), sg.Combo(available_taxonomic_levels_list, default_value="Species", key="parcat_taxonomic_level")]], title="Settings")],

    					[sg.Text('',size=(1,1))],
    					[sg.Text('Species occurrence analysis',size=(40,2), font=('Arial', 11, "bold"))],
                        [sg.Text('Collect species occurrences from GBIF'), sg.Button('Download', key="run_species_occurrence"), sg.Button('Plot', key="run_species_occurrence_plot"),
                        sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("1000", size=(4,1), key="width_occurrence"),
                        sg.Input("1000", size=(4,1), key="height_occurrence")]], title="Settings")],
                        [sg.Frame(layout=[[sg.CB("Africa", key="africa_occurrence"), sg.CB("Antarctica", key="antarctica_occurrence"),
                        sg.CB("Europe", key="europe_occurrence"), sg.CB("Asia", key="asia_occurrence"),
                        sg.CB("North America", key="north_america_occurrence"), sg.CB("Oceania", key="oceania_occurrence"),
                        sg.CB("South America", key="south_america_occurrence")]], title="Continents to check")],
    					]

    read_proportions_layout =  [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Read proportions', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(initial_folder = taXon_table_dir_path, key = 'taxon_table_7_path'), sg.Button("Open", key="open7")],
    					[sg.Text('_'*120)],
    					[sg.Text('',size=(1,1))],

    					[sg.Text('Read proportions',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Calculate heatmap", size=(25,1)), sg.Button("Run", key = 'run_read_proportions_heatmap'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h): "), sg.Input("1000", size=(4,1), key="width_read_props_heatmap"),
                        sg.Input("1000", size=(4,1), key="heigth_read_props_heatmap"), sg.Text("Tax. level:"), sg.Combo(available_taxonomic_levels_list, default_value="Class", key="read_proportions_heatmap_taxonomic_level")]], title="Settings")],
                        [sg.Text("Calculate bar chart", size=(25,1)), sg.Button("Run", key = 'run_read_proportions_bar'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h): "), sg.Input("1000", size=(4,1), key="width_read_props_bar"),
                        sg.Input("1000", size=(4,1), key="heigth_read_props_bar"), sg.Text("Tax. level:"), sg.Combo(available_taxonomic_levels_list, default_value="Class", key="read_proportions_bar_taxonomic_level")]], title="Settings"), sg.Button(" ? ", key="run_read_proportions_help", button_color=("Blue", "White"))],
                        [sg.Text("Calculate pie chart", size=(25,1)), sg.Button("Run", key = 'run_read_proportions_pie'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h): "), sg.Input("1000", size=(4,1), key="width_read_props_pie"),
                        sg.Input("1000", size=(4,1), key="heigth_read_props_pie"), sg.Text("Tax. level:"), sg.Combo(available_taxonomic_levels_list, default_value="Class", key="read_proportions_pie_taxonomic_level")]], title="Settings")],
    					[sg.Text('',size=(1,1))],

    					[sg.Text('Krona charts',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Create Krona chart", size=(25,1)), sg.Button("Run", key = 'run_create_krona_chart'), sg.Text("", size=(1,1)),
                        sg.Frame(layout=[[sg.Text("Layout:"), sg.Combo(["Merged data", "Single samples"], default_value="Single samples", key="krona_version")],
                        [sg.Text("Charts are based on read proportions. Krona tools is required.")]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					]

    ## load available colorscales
    heatmap_colorscale_list = px.colors.named_colorscales()

    diversity_analysis_layout =  [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Diversity analyses', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(initial_folder = taXon_table_dir_path, key = 'taxon_table_8_path'), sg.Button("Open", key="open8")],
    					[sg.Text('_'*120)],
    					[sg.Text('',size=(1,1))],

    					[sg.Text('Alpha diversity',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Calculate observed taxa:", size=(28,1)), sg.Button("Run", key = 'run_alpha_diversity'),
                        sg.Button(" ? ", key = 'run_alpha_diversity_help_text', button_color=('blue', 'white')),
                        sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):", size=(11,1)), sg.Input("1000", size=(4,1), key="width_alpha"),
                        sg.Input("1000", size=(4,1), key="height_alpha"), sg.Text('Scatter size:'), sg.Input("15", size=(4,1), key="alpha_scatter_size")],
                        [sg.Text('Plot type:'), sg.Combo(["Scatter plot", "Box plot"], default_value="Box plot", key="alpha_diversity_plot"),
                        sg.Text('Tax. level:'), sg.Combo(available_taxonomic_levels_list, default_value='OTUs', key="alpha_diversity_taxonomic_level")]], title="Settings")],
    					[sg.Text('',size=(1,1))],

    					[sg.Text('Beta diversity',size=(40,2), font=('Arial', 11, "bold"))],
                        [sg.Text("Calculate beta diversity cluster:", size=(28,1)), sg.Button("Run", key = 'run_betadiv_clustering'), sg.Button(" ? ", key = 'run_cluster_dendrogram_help_text', button_color=('blue', 'white')), sg.Text("", size=(1,1)),
                        sg.Frame(layout=[[sg.Text("Plot size (w,h):", size=(11,1)), sg.Input("1000", size=(4,1), key="width_betadiv_clustering"),
                        sg.Input("1000", size=(4,1), key="height_betadiv_clustering"), sg.Text('Tax. level:'), sg.Combo(available_taxonomic_levels_list, default_value='OTUs', key="betadiv_cluster_taxonomic_level")],
                        [sg.Text('Linkage:', size=(11,1)), sg.Combo(['single', 'average', 'complete'], default_value='average', size=(8,1), key="betadiv_linkage"), sg.Text('Group threshold:'), sg.Combo([str(round(i,1)) for i in np.arange(0.1, 1.0, 0.1)], default_value='0.3', size=(5,1), key="betadiv_clustering_threshold")],
                        [sg.Text("Method:"), sg.Combo(["jaccard", "braycurtis"], default_value="jaccard", key="betadiv_cluster_dissimilarity")]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text("Calculate beta diversity heatmap:", size=(28,1)), sg.Button("Run", key = 'run_beta_diversity'),
                        sg.Button(" ? ", key = 'run_beta_diversity_help_text', button_color=('blue', 'white')),
                        sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):", size=(11,1)), sg.Input("1000", size=(4,1), key="width_beta"),
                        sg.Input("1000", size=(4,1), key="height_beta"), sg.Text('Tax. level:'), sg.Combo(available_taxonomic_levels_list, default_value='OTUs', key="betadiv_heatmap_taxonomic_level")],
                        [sg.Text("Method:"), sg.Combo(["jaccard", "braycurtis"], default_value="jaccard", key="betadiv_heatmap_dissimilarity")]], title="Settings")],
    					[sg.Text('',size=(1,1))],
        				]

    ordination_analysis_layout =  [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Ordination analyses', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(initial_folder = taXon_table_dir_path, key = 'taxon_table_9_path'), sg.Button("Open", key="open9")],
    					[sg.Text('_'*120)],
    					[sg.Text('',size=(1,1))],

    					[sg.Text('Principle coordinate analysis',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Calculate PCoA", size=(23,1)), sg.Button("Run", key = 'run_PCoA_analysis'),
                        sg.Button(" ? ", key = 'run_pcoa_analysis_help_text', button_color=('blue', 'white')),
                        sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("1000", size=(4,1), key="width_pcoa"),
                        sg.Input("1000", size=(4,1), key="height_pcoa"), sg.Text("Scatter size:"), sg.Input("10", size=(3,1), key="scatter_size_pcoa"),
                        sg.Text('Tax. level:'), sg.Combo(available_taxonomic_levels_list, default_value='OTUs', key="pcoa_analysis_taxonomic_level")],
                        [sg.Text("Method:"), sg.Combo(["jaccard", "braycurtis"], default_value="jaccard", key="pcoa_dissimilarity")]], title="Settings")],
    					[sg.Text('',size=(1,1))],

    					[sg.Text('Non-metric multidimensional scaling analysis',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Calculate NMDS", size=(23,1)), sg.Button("Run", key = 'run_NMDS_analysis'),
                        sg.Button(" ? ", key = 'run_nmds_analysis_help_text', button_color=('blue', 'white')),
                        sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("1000", size=(4,1), key="width_nmds"),
                        sg.Input("1000", size=(4,1), key="height_nmds"), sg.Text("Scatter size:"), sg.Input("10", size=(3,1), key="nmds_s"),
                        sg.Text('Tax. level:'), sg.Combo(available_taxonomic_levels_list, default_value='OTUs', key="nmds_analysis_taxonomic_level")],
                        [sg.Text("Method:"), sg.Combo(["jaccard", "braycurtis"], default_value="jaccard", key="nmds_dissimilarity"),
                        sg.Text("Max. iterations:"), sg.Input("1000", size=(4,1), key="max_iter_val"), sg.Text("No. of initiations:"), sg.Input("100", size=(4,1), key="n_init_val")]], title="Settings")],
    					[sg.Text('',size=(1,1))],

    					[sg.Text('Canonical coordination analysis',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Calculate CCA", size=(23,1)), sg.Button("Run", key = 'run_CCA_analysis'),
                        sg.Button(" ? ", key = 'run_cca_analysis_help_text', button_color=('blue', 'white')),
                        sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("1000", size=(4,1), key="width_cca"),
                        sg.Input("1000", size=(4,1), key="height_cca"), sg.Text("Scatter size:"), sg.Input("10", size=(3,1), key="cca_scatter_size"),
                        sg.Text('Tax. level:'), sg.Combo(available_taxonomic_levels_list, default_value='OTUs', key="cca_analysis_taxonomic_level")]], title="Settings")],
    					[sg.Text('',size=(1,1))],

        				]

    taxon_list_layout = [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Taxon list',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(initial_folder = taXon_table_dir_path, key = 'taxon_table_10_path'), sg.Button("Open", key="open10")],
    					[sg.Text('_'*120)],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Convert to taxon list',size=(40,2), font=('Arial', 11, "bold"))],
                        [sg.Text("Enter name of output file:"), sg.Input(project_folder + "_taxon_list", size=(48,1), key='taxon_list_output_file_name'), sg.Text(".txt + .xlsx")],
                        [sg.Text('GBIF link:', size=(20,1)), sg.CB("", key="create_gbif_link", default=True), sg.Frame(layout=[[sg.Text("Create link to species entries and search for vernacular names in the GBIF database.")]], title="Information")],
                        [sg.Text('Intraspecific distances:', size=(20,1)), sg.CB("", key="calc_dist", default=True), sg.Frame(layout=[[sg.Text("Calculate intraspecific distances for species represented by more than one OTU.")]], title="Information")],
                        [sg.Text('Use metadata:', size=(20,1)), sg.CB("", key="use_metadata", default=True), sg.Frame(layout=[[sg.Text("Calculate the taxon occurruence according to the metadata.")]], title="Information")],
                        [sg.Text('', size=(1,1))],
                        [sg.Text("Create taxon list"), sg.Button("Run", key = 'run_create_taxon_list')],
    					[sg.Text('',size=(1,1))]
                        ]

    create_perlodes_input_layout = [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Water Framework Directive (WFD)',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(initial_folder = taXon_table_dir_path, key = 'taxon_table_11_path'), sg.Button("Open", key="open11")],
    					[sg.Text('_'*120)],
    					[sg.Text('',size=(1,1))],
        				[sg.Text('Create Perlodes input file',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text('Operational taxon list:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'operational_taxon_list_path'),
                        sg.Button('Download', key='open_operational_taxonlist_download', button_color=('black', 'white'))],
    					[sg.Text('',size=(1,1))],
    					[sg.Text("Create Perlodes input file: "), sg.Button("Run", key = 'run_create_perlodes_input_file')],
    					[sg.Text('',size=(1,2))],
    					[sg.Button(key = 'open_fgbewertung', button_color=('white', 'white'), image_filename=fliessgewaesserbewertungde)],
                        ]

    layout = [  [sg.Image(ttt_logo), sg.Text("", size=(9,1)), sg.Text('Project:', font=('Arial', 12, "bold")), sg.Text(project_folder, font=('Arial', 12, "bold"))],
    			[sg.Text('',size=(1,1))],
    			[sg.TabGroup([[
                sg.Tab('Start\n', getting_started_layout),
                sg.Tab('Replicate\nanalyses', replicate_analyses_layout),
                sg.Tab('Table\nfiltering', taxon_table_filtering_layout),
                sg.Tab('Data\nconversion', data_conversion_layout),
                sg.Tab('Basic\nstatistics', basic_statistics_layout),
                sg.Tab('Sample\ncomparison', sample_comparison_layout),
                sg.Tab('Taxonomic\ncomposition', taxonomic_composition_layout),
                sg.Tab('Read\nproportions', read_proportions_layout),
                sg.Tab('Diversity\nanalyses', diversity_analysis_layout),
                sg.Tab('Ordination\nanalyses', ordination_analysis_layout),
                sg.Tab('Taxon\nlist', taxon_list_layout),
                sg.Tab('Water Framework\nDirective', create_perlodes_input_layout)]])],
    			[sg.Text('',size=(1,1))],
    			[sg.Exit(button_color=('black', 'red')), sg.Text("", size=(35,1)),
                sg.Button("Open log file", button_color=('black', 'white'), key="run_open_log_file"), sg.Text("", size=(3,1)),
                sg.Button("Log file network", button_color=('black', 'white'), key='run_log_network'), sg.Text("", size=(3,1)),
                sg.Button("Tutorial", button_color=('black', 'white'), key="open_tutorial"), sg.Text("", size=(3,1)),
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
            ## collect discete color sequence
            color_discrete_sequence = available_colorsequences[values["colorsequence"]]
            ## collect the color scale
            colorscale = values["colorscale"]
            ## collect plotly template
            template = values["template"]
            ## collect opactity value (transparent theme)
            if values['transparent_theme'] != True:
                opacity_value = 1
            else:
                opacity_value = 0.5
            theme = [values["color1"], values["color2"], opacity_value]
            ## collect font size
            font_size = values["font_size"]

            # define variables
            TaXon_table_name = values['TaXon_table_name']
            taxon_table_1_path = values['taxon_table_1_path']
            replicate_suffix_entry = values['replicate_suffix_entry']
            taxonomy_table_path = values['taxonomy_table_path']
            Read_table_path = values['Read_table_path']
            venn_taxon_table_2_path = values['venn_taxon_table_2_path']
            venn_taxon_table_3_path = values['venn_taxon_table_3_path']
            venn_diagram_name = values['venn_diagram_name']
            taxon_table_8_path = values['taxon_table_8_path']
            taxon_table_2_path = values['taxon_table_2_path']
            taxon_table_3_path = values['taxon_table_3_path']
            filter_taxonomic_level = values['filter_taxonomic_level']
            taxon_table_4_path = values['taxon_table_4_path']
            taxon_table_5_path = values['taxon_table_5_path']
            taxon_table_6_path = values['taxon_table_6_path']
            taxon_table_7_path = values['taxon_table_7_path']
            taxon_list_output_file_name = values['taxon_list_output_file_name']
            sample_repetitions = values['sample_repetitions']
            read_repetitions = values['read_repetitions']
            width_tax_res = values["width_tax_res"]
            height_tax_res = values["height_tax_res"]
            width_tax_rich = values["width_tax_rich"]
            height_tax_rich = values["height_tax_rich"]
            heigth_read_props_heatmap = values["heigth_read_props_heatmap"]
            heigth_read_props_bar = values["heigth_read_props_bar"]
            heigth_read_props_pie = values["heigth_read_props_pie"]
            width_read_props_heatmap = values["width_read_props_heatmap"]
            width_read_props_bar = values["width_read_props_bar"]
            width_read_props_pie = values["width_read_props_pie"]
            width_site_occ = values["width_site_occ"]
            height_site_occ = values["height_site_occ"]
            operational_taxon_list_path = values["operational_taxon_list_path"]
            taxon_table_9_path = values["taxon_table_9_path"]
            taxon_table_10_path = values["taxon_table_10_path"]
            taxon_table_11_path = values["taxon_table_11_path"]
            width_alpha = values["width_alpha"]
            height_alpha = values["height_alpha"]
            alpha_scatter_size = values["alpha_scatter_size"]
            width_beta = values["width_beta"]
            height_beta = values["height_beta"]
            tres_type_a = values["tres_type_a"]
            tres_type_b = values["tres_type_b"]
            cm_taxon_based_filter_all = values['cm_taxon_based_filter_all']
            cm_taxon_based_filter_none = values['cm_taxon_based_filter_none']
            cm_sample_based_filter_all = values['cm_sample_based_filter_all']
            cm_sample_based_filter_none = values['cm_sample_based_filter_none']
            taxon_filter_method_keep = values['taxon_filter_method_keep']
            taxon_filter_method_exclude = values['taxon_filter_method_exclude']
            sample_filter_method_keep = values['sample_filter_method_keep']
            sample_filter_method_exclude = values['sample_filter_method_exclude']
            sample_filter_method_criterion = values['sample_filter_method_criterion']
            taxon_filter_method_criterion = values['taxon_filter_method_criterion']
            create_gbif_link = values['create_gbif_link']
            calc_dist = values['calc_dist']
            use_metadata = values['use_metadata']
            sorting_method_fh = values['sorting_method_fh']
            sorting_method_jamp = values['sorting_method_jamp']
            sorting_method_boldigger = values['sorting_method_boldigger']
            width_cca = values['width_cca']
            height_cca = values['height_cca']
            cca_scatter_size = values['cca_scatter_size']
            width_pcoa = values['width_pcoa']
            height_pcoa = values['height_pcoa']
            scatter_size_pcoa = values['scatter_size_pcoa']
            height_bstats = values["height_bstats"]
            width_bstats = values["width_bstats"]
            taxon_table_0_path = values["taxon_table_0_path"]
            read_table_format_TTT = values["read_table_format_TTT"]
            read_table_format_qiime = values["read_table_format_qiime"]
            height_repicates = values["height_repicates"]
            width_repicates = values["width_repicates"]
            width_repcorr = values["width_repcorr"]
            height_repcorr = values["height_repcorr"]
            height_per_taxon_stats = values["height_per_taxon_stats"]
            width_per_taxon_stats = values["width_per_taxon_stats"]
            per_taxon_stats_taxonomic_level = values["per_taxon_stats_taxonomic_level"]
            absolute_filtering = values["absolute_filtering"]
            absolute_filter_threshold = values["absolute_filter_threshold"]
            relative_filtering = values["relative_filtering"]
            relative_filter_threshold = values["relative_filter_threshold"]
            height_rarefaction_reads = values["height_rarefaction_reads"]
            width_rarefaction_reads = values["width_rarefaction_reads"]
            rarefaction_curve_taxonomic_level_1 = values["rarefaction_curve_taxonomic_level_1"]
            rarefaction_curve_taxonomic_level_2 = values["rarefaction_curve_taxonomic_level_2"]
            y_zero_repcorr = values["y_zero_repcorr"]
            x_zero_repcorr = values["x_zero_repcorr"]
            africa_occurrence = values["africa_occurrence"]
            antarctica_occurrence = values["antarctica_occurrence"]
            asia_occurrence = values["asia_occurrence"]
            europe_occurrence = values["europe_occurrence"]
            north_america_occurrence = values["north_america_occurrence"]
            oceania_occurrence = values["oceania_occurrence"]
            south_america_occurrence = values["south_america_occurrence"]
            height_occurrence = values["height_occurrence"]
            width_occurrence = values["width_occurrence"]
            height_betadiv_clustering = values["height_betadiv_clustering"]
            width_betadiv_clustering = values["width_betadiv_clustering"]
            betadiv_clustering_threshold = values["betadiv_clustering_threshold"]
            betadiv_linkage = values["betadiv_linkage"]
            site_occupancy_taxonomic_level = values["site_occupancy_taxonomic_level"]
            read_proportions_heatmap_taxonomic_level = values["read_proportions_heatmap_taxonomic_level"]
            read_proportions_bar_taxonomic_level = values["read_proportions_bar_taxonomic_level"]
            read_proportions_pie_taxonomic_level = values["read_proportions_pie_taxonomic_level"]
            betadiv_heatmap_taxonomic_level = values["betadiv_heatmap_taxonomic_level"]
            betadiv_cluster_taxonomic_level = values["betadiv_cluster_taxonomic_level"]
            alpha_diversity_taxonomic_level = values["alpha_diversity_taxonomic_level"]
            alpha_diversity_plot = values["alpha_diversity_plot"]
            pcoa_analysis_taxonomic_level = values["pcoa_analysis_taxonomic_level"]
            cca_analysis_taxonomic_level = values["cca_analysis_taxonomic_level"]
            nmds_analysis_taxonomic_level = values["nmds_analysis_taxonomic_level"]
            parcat_taxonomic_level = values["parcat_taxonomic_level"]
            width_parcat = values["width_parcat"]
            height_parcat = values["height_parcat"]
            nmds_s = values["nmds_s"]
            max_iter_val = values["max_iter_val"]
            n_init_val = values["n_init_val"]
            width_nmds = values["width_nmds"]
            height_nmds = values["height_nmds"]
            cm_nc_subtraction_all = values["cm_nc_subtraction_all"]
            nc_subtraction_criterion = values["nc_subtraction_criterion"]
            consistency_option = values["consistency_option"]
            nmds_dissimilarity = values["nmds_dissimilarity"]
            pcoa_dissimilarity = values["pcoa_dissimilarity"]
            betadiv_cluster_dissimilarity = values["betadiv_cluster_dissimilarity"]
            betadiv_heatmap_dissimilarity = values["betadiv_heatmap_dissimilarity"]
            rarefaction_curve_method = values["rarefaction_curve_method"]
            krona_version = values["krona_version"]
            replicate_analysis_custom_colors = values["replicate_analysis_custom_colors"]
            clustering_unit = values["clustering_unit"]
            site_occupancy_option = values["site_occupancy_option"]
            site_occupancy_heatmap_add_metadata_sum = values["site_occupancy_heatmap_add_metadata_sum"]

            if event == 'run_taxon_converter':

                if (taxonomy_table_path == '' or Read_table_path == ''):
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    if sorting_method_fh == True:
                        sheet_name = "First hit"
                    elif sorting_method_jamp == True:
                        sheet_name = "JAMP hit"
                    else:
                        sheet_name = "BOLDigger hit"

                    if read_table_format_TTT == True:
                        from taxontabletools.taxon_table_converter import taxon_table_converter_ttt
                        taxon_table_converter_ttt(Read_table_path, taxonomy_table_path, TaXon_table_name, sheet_name, path_to_outdirs)
                    elif read_table_format_qiime == True:
                        from taxontabletools.taxon_table_converter import taxon_table_converter_qiime2
                        taxon_table_converter_qiime2(Read_table_path, taxonomy_table_path, TaXon_table_name, sheet_name, path_to_outdirs)

            if event == 'run_taxonomy_table_check':
                if (taxonomy_table_path == ''):
                    sg.PopupError("Please provide a file", keep_on_top=True)
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
                if (Read_table_path == ''):
                    sg.PopupError("Please provide a file", keep_on_top=True)
                elif read_table_format_TTT == True:
                    from taxontabletools.check_read_table_format import check_read_table_format_TTT
                    check_read_table_format_TTT(Read_table_path)
                elif read_table_format_qiime == True:
                    from taxontabletools.check_read_table_format import check_read_table_format_qiime2
                    check_read_table_format_qiime2(Read_table_path)

            if event == 'run_taXon_table_check':
                if (taxon_table_0_path == ''):
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.check_taXon_table_format import check_taXon_table_format
                    check_taXon_table_format(taxon_table_0_path)

            if event == "run_taxon_table_per_sample":
                if taxon_table_3_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.create_taxon_table_per_sample import create_taxon_table_per_sample
                    create_taxon_table_per_sample(taxon_table_3_path, path_to_outdirs)

            if event == "run_convert_to_pa":
                if taxon_table_3_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.convert_to_presence_absence import convert_to_presence_absence
                    convert_to_presence_absence(taxon_table_3_path, path_to_outdirs)

            if event == 'run_replicate_analysis':
                if taxon_table_1_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    replicate_suffix_list = replicate_suffix_entry.replace(" ", "").split(",")
                    if replicate_analysis_custom_colors == True:
                        layout2 = [[sg.Text("Custom colors for plot 3:", size=(20,1))],
                        [sg.Text('Shared OTUs:', size=(15,1)), sg.Combo(plotly_colors, key="c1", default_value="Teal")],
                        [sg.Text("Non shared OTUs:", size=(15,1)), sg.Combo(plotly_colors, key="c2", default_value="Orange")],
                        [sg.Button('Ok')],
                        [sg.Text("")]]
                        win2 = sg.Window('Replicate analysis options', layout2, grab_anywhere=True)
                        while True:
                            event2, values2 = win2.Read()
                            custom_colors = [values2["c1"], values2["c2"]]
                            win2.Close()
                            win2_active = False
                            break
                    else:
                        custom_colors = color_discrete_sequence
                    from taxontabletools.replicate_analysis import replicate_analysis
                    replicate_analysis(taxon_table_1_path, height_repicates, width_repicates, replicate_suffix_list, path_to_outdirs, template, theme, font_size, custom_colors, clustering_unit)

            if event == 'run_replicate_correlation':
                if taxon_table_1_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    repcorr_options = [x_zero_repcorr, y_zero_repcorr]
                    replicate_suffix_list = replicate_suffix_entry.replace(" ", "").split(",")
                    from taxontabletools.replicate_correlation_analysis import replicate_correlation_analysis
                    replicate_correlation_analysis(taxon_table_1_path, replicate_suffix_list, path_to_outdirs, width_repcorr, height_repcorr, repcorr_options, template, theme, font_size, clustering_unit)

            if event == 'run_venn2_diagram':
                if os.path.exists(str(path_to_outdirs) + "/Venn_diagrams/" + venn_diagram_name):
                    sg.PopupError("Error: Venn directory already exists", keep_on_top=True)
                elif (taxon_table_5_path == '' or venn_taxon_table_2_path == ''):
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.venn_diagram import venn_diagram
                    venn_taxon_table_3_path = False
                    venn_diagram(taxon_table_5_path, venn_taxon_table_2_path, venn_taxon_table_3_path, venn_diagram_name, path_to_outdirs, clustering_unit)

            if event == 'run_venn3_diagram':
                if os.path.exists(str(path_to_outdirs) + "/Venn_diagrams/" + venn_diagram_name):
                    sg.PopupError("Error: Venn directory already exists", keep_on_top=True)
                elif (taxon_table_5_path == '' or venn_taxon_table_2_path == ''):
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.venn_diagram import venn_diagram
                    venn_diagram(taxon_table_5_path, venn_taxon_table_2_path, venn_taxon_table_3_path, venn_diagram_name, path_to_outdirs, clustering_unit)

            if event == 'run_rarefaction_curve':
                if taxon_table_5_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    if rarefaction_curve_method == "all-in-one":
                        from taxontabletools.rarefaction_curve import rarefaction_curve_legacy
                        if rarefaction_curve_taxonomic_level_1 == "OTUs":
                            rarefaction_curve_taxonomic_level_1 = clustering_unit
                        rarefaction_curve_legacy(taxon_table_5_path, int(sample_repetitions), path_to_outdirs, template, theme, font_size, rarefaction_curve_taxonomic_level_1)
                    else:
                        if rarefaction_curve_taxonomic_level_1 == "OTUs":
                            rarefaction_curve_taxonomic_level_1 = clustering_unit
                        if rarefaction_curve_taxonomic_level_2 == "OTUs":
                            rarefaction_curve_taxonomic_level_2 = clustering_unit
                        from taxontabletools.rarefaction_curve import rarefaction_curve_taxa
                        rarefaction_curve_taxa(taxon_table_5_path, int(sample_repetitions), path_to_outdirs, template, font_size, rarefaction_curve_taxonomic_level_1, rarefaction_curve_taxonomic_level_2, color_discrete_sequence)

            if event == 'run_rarefaction_curve_reads':
                if taxon_table_5_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    answer = sg.PopupOKCancel("This analysis can take very long. Continue?", title="Runtime")
                    if answer == "OK":
                        from taxontabletools.rarefaction_curve import rarefaction_curve_reads
                        rarefaction_curve_reads(taxon_table_5_path, int(read_repetitions), width_rarefaction_reads, height_rarefaction_reads, path_to_outdirs, template, theme, font_size)

            if event == 'run_create_taxon_list':
                if taxon_table_10_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.create_taxon_list import create_taxon_list
                    create_taxon_list(taxon_table_10_path, taxon_list_output_file_name, create_gbif_link, calc_dist, use_metadata, taxon_tools_version, path_to_outdirs, clustering_unit)

            if event == 'run_taxon_filtering' and not win2_active:

                if taxon_table_2_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.taxon_table_filtering import create_taxa_mask
                    available_taxa = create_taxa_mask(taxon_table_2_path, filter_taxonomic_level, taxon_filter_method_criterion)
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
                    [sg.Frame(layout = available_taxa_list, title = 'Check taxa to filter')],
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

                            from taxontabletools.taxon_table_filtering import taxon_filter
                            if taxon_filter_method_keep == True:
                                taxon_filter_method = "keep"
                            else:
                                taxon_filter_method = "exclude"
                            if check_list == []:
                                check_list = ['no taxa to filter']
                                taxon_filter(taxon_table_2_path, check_list, filter_taxonomic_level, appendix_name, threshold, path_to_outdirs, taxon_filter_method)
                            else:
                                taxon_filter(taxon_table_2_path, check_list, filter_taxonomic_level, appendix_name, threshold, path_to_outdirs, taxon_filter_method)

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
                if taxon_table_2_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.taxon_table_filtering import create_sample_mask
                    available_samples = sorted(create_sample_mask(taxon_table_2_path, sample_filter_method_criterion))
                    win2_active = True
                    window.Hide()
                    if cm_sample_based_filter_all == True:
                        available_samples_list = list(slices([sg.CB(name, default=True) for name in sorted(available_samples)], 10))
                    else:
                        available_samples_list = list(slices([sg.CB(name) for name in sorted(available_samples)], 10))
                    layout2 = [[sg.Text("Sample filter", size=(20,1))],
                    [sg.Frame(layout = available_samples_list, title = 'Check samples to filter')],
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

                            from taxontabletools.taxon_table_filtering import filter_samples
                            if sample_filter_method_keep == True:
                                sample_filter_method = "keep"
                            else:
                                sample_filter_method = "exclude"
                            filter_samples(taxon_table_2_path, check_list, appendix_name, path_to_outdirs, sample_filter_method)

                            win2.Close()
                            win2_active = False
                            window.UnHide()
                            break

                        if event2 is None or event2 == 'Back':
                            win2.Close()
                            win2_active = False
                            window.UnHide()
                            break

            if event == 'run_nc_subtraction' and not win2_active:
                if taxon_table_2_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.taxon_table_filtering import create_sample_mask
                    available_samples = sorted(create_sample_mask(taxon_table_2_path, nc_subtraction_criterion))
                    win2_active = True
                    window.Hide()
                    if cm_nc_subtraction_all == True:
                        available_samples_list = list(slices([sg.CB(name, default=True) for name in sorted(available_samples)], 10))
                    else:
                        available_samples_list = list(slices([sg.CB(name) for name in sorted(available_samples)], 10))
                    layout2 = [[sg.Text("Subtraction of negative controls", size=(20,1))],
                    [sg.Frame(layout = available_samples_list, title = 'Select your negative controls')],
                    [sg.Button('Subtract')],
                    [sg.Button('Back')]]

                    win2 = sg.Window('TaXon table filtering', layout2, keep_on_top=True, grab_anywhere=True)

                    while True:
                        event2, values2 = win2.Read()

                        if event2 == 'Subtract':
                            i = 0
                            check_dict = {}
                            check_list = []

                            for key, value in values2.items():
                                check_dict[available_samples[i]] = value
                                i += 1
                            for key, value in check_dict.items():
                                if value == True:
                                    check_list.append(key)

                            from taxontabletools.taxon_table_filtering import subtract_NCs
                            subtract_NCs(taxon_table_2_path, path_to_outdirs, check_list)

                            win2.Close()
                            win2_active = False
                            window.UnHide()
                            break

                        if event2 is None or event2 == 'Back':
                            win2.Close()
                            win2_active = False
                            window.UnHide()
                            break

            if event == 'run_gbif_check_taxonomy':
                if taxon_table_3_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.gbif_check_taxonomy import gbif_check_taxonomy
                    gbif_check_taxonomy(taxon_table_3_path, path_to_outdirs)

            if event == 'run_preview_suffix':
                replicate_suffix_list = replicate_suffix_entry.replace(" ", "").split(",")
                text = "\n" + "Sample_No1_comb = " + ', '.join(["Sample_No1_" + suffix for suffix in replicate_suffix_list]) + "\n"
                sg.Popup(text, title='Suffix preview')

            if event == 'run_taxonomic_richness':
                if taxon_table_4_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.calculate_taxonomic_richness import calculate_taxonomic_richness
                    calculate_taxonomic_richness(taxon_table_4_path, path_to_outdirs, width_tax_rich, height_tax_rich, template, theme, font_size)

            if event == 'run_taxonomic_resolution':
                if taxon_table_4_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    if tres_type_a == True:
                        figure_type = "a"
                    if tres_type_b == True:
                        figure_type = "b"
                    from taxontabletools.calculate_taxonomic_resolution import calculate_taxonomic_resolution
                    calculate_taxonomic_resolution(taxon_table_4_path, path_to_outdirs, width_tax_res, height_tax_res, figure_type, template, theme, font_size, clustering_unit)

            if event == 'run_create_meta_data_table':
                if taxon_table_3_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.create_metadata_table import create_metadata_table
                    create_metadata_table(taxon_table_3_path, path_to_outdirs)

            if event == "run_modify_metadata_table":
                    from taxontabletools.create_metadata_table import modify_metadata_table
                    modify_metadata_table(path_to_outdirs)

            if event == 'run_site_occupancy':
                if taxon_table_5_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:

                    if site_occupancy_taxonomic_level == "OTUs":
                        site_occupancy_taxonomic_level = clustering_unit

                    from taxontabletools.get_available_meta_data import get_available_meta_data
                    meta_data_to_test = get_available_meta_data(taxon_table_5_path, path_to_outdirs)

                    if meta_data_to_test == False:
                                 sg.PopupError("Metadata table not found for the loaded TaXon table!\nPlease first create a metadata table!", keep_on_top=True)
                    else:
                                 win2_active = True
                                 window.Hide()
                                 meta_data_to_test_list = list(slices([sg.Radio(name, "metadata", key=name, default=True) for name in sorted(meta_data_to_test)], 5))
                                 layout2 = [[sg.Text("Site occupancy", size=(20,1))],
                                 [sg.Frame(layout = meta_data_to_test_list, title = 'Check metadata to test')],
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
                                        if site_occupancy_option == "barchart":
                                            from taxontabletools.site_occupancy import site_occupancy_barchart
                                            site_occupancy_barchart(taxon_table_5_path, meta_data_to_test, site_occupancy_taxonomic_level, path_to_outdirs, width_site_occ, height_site_occ, template, theme, font_size)
                                            win2.Close()
                                            win2_active = False
                                            window.UnHide()
                                        else:
                                            from taxontabletools.site_occupancy import site_occupancy_heatmap
                                            site_occupancy_heatmap(taxon_table_5_path, path_to_outdirs, template, height_site_occ, width_site_occ, meta_data_to_test, site_occupancy_taxonomic_level, font_size, color_discrete_sequence, site_occupancy_heatmap_add_metadata_sum)
                                            win2.Close()
                                            win2_active = False
                                            window.UnHide()

                                    if event2 is None or event2 == 'Back':
                                        win2.Close()
                                        win2_active = False
                                        window.UnHide()
                                        break

            if event == 'run_read_proportions_heatmap':
                if taxon_table_7_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    if read_proportions_heatmap_taxonomic_level == "OTUs":
                        read_proportions_heatmap_taxonomic_level = clustering_unit
                    from taxontabletools.read_proportions import read_proportions_heatmap
                    read_proportions_heatmap(taxon_table_7_path, read_proportions_heatmap_taxonomic_level, path_to_outdirs, width_read_props_heatmap, heigth_read_props_heatmap, template, font_size)

            if event == 'run_read_proportions_bar':
                if taxon_table_7_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    if read_proportions_bar_taxonomic_level == "OTUs":
                        read_proportions_bar_taxonomic_level = clustering_unit
                    from taxontabletools.read_proportions import read_proportions_bar
                    read_proportions_bar(taxon_table_7_path, read_proportions_bar_taxonomic_level, path_to_outdirs, width_read_props_bar, heigth_read_props_bar, template, font_size, color_discrete_sequence)

            if event == 'run_read_proportions_pie':
                if taxon_table_7_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    if read_proportions_pie_taxonomic_level == "OTUs":
                        read_proportions_pie_taxonomic_level = clustering_unit
                    from taxontabletools.read_proportions import read_proportions_pie
                    read_proportions_pie(taxon_table_7_path, read_proportions_pie_taxonomic_level, path_to_outdirs, width_read_props_pie, heigth_read_props_pie, template, font_size, color_discrete_sequence)

            if event == 'run_basic_stats':
                if taxon_table_4_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.basic_stats import basic_stats
                    basic_stats(taxon_table_4_path, height_bstats, width_bstats, path_to_outdirs, template, theme, font_size, clustering_unit)

            if event == 'run_create_krona_chart':
                if taxon_table_7_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                elif krona_version == "Merged data":
                    from taxontabletools.create_krona_chart import create_krona_chart_single
                    create_krona_chart_single(taxon_table_7_path, path_to_outdirs)
                else:
                    from taxontabletools.create_krona_chart import create_krona_chart_multi
                    create_krona_chart_multi(taxon_table_7_path, path_to_outdirs)

            if event == 'run_CCA_analysis':
                if taxon_table_9_path == '':
                    sg.PopupError("Please provide a file")
                else:

                    if cca_analysis_taxonomic_level == "OTUs":
                        cca_analysis_taxonomic_level = clustering_unit

                    from taxontabletools.get_available_meta_data import get_available_meta_data
                    meta_data_to_test = get_available_meta_data(taxon_table_9_path, path_to_outdirs)

                    if meta_data_to_test == False:
                            sg.PopupError("Metadata table not found for the loaded TaXon table!\nPlease first create a metadata table!")
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
                                                        CCA_analysis(taxon_table_9_path, meta_data_to_test, cca_analysis_taxonomic_level, width_cca, height_cca, cca_scatter_size, path_to_outdirs, template, font_size, color_discrete_sequence)
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

                if taxon_table_9_path == '':
                    sg.PopupError("Please provide a file")
                else:

                    if pcoa_analysis_taxonomic_level == "OTUs":
                        pcoa_analysis_taxonomic_level = clustering_unit

                    from taxontabletools.get_available_meta_data import get_available_meta_data
                    meta_data_to_test = get_available_meta_data(taxon_table_9_path, path_to_outdirs)

                    if meta_data_to_test == False:
                            sg.PopupError("Metadata table not found for the loaded TaXon table!\nPlease first create a metadata table!")
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
                                                        PCoA_analysis(taxon_table_9_path, meta_data_to_test, pcoa_analysis_taxonomic_level, width_pcoa, height_pcoa, scatter_size_pcoa, path_to_outdirs, template, font_size, color_discrete_sequence, pcoa_dissimilarity)
                                                        win2.Close()
                                                        win2_active = False
                                                        window.UnHide()
                                                        break

                                if event2 is None or event2 == 'Back':
                                        win2.Close()
                                        win2_active = False
                                        window.UnHide()
                                        break

            if event == 'run_NMDS_analysis':

                if taxon_table_9_path == '':
                    sg.PopupError("Please provide a file")
                else:

                    if nmds_analysis_taxonomic_level == "OTUs":
                        nmds_analysis_taxonomic_level = clustering_unit

                    from taxontabletools.get_available_meta_data import get_available_meta_data
                    meta_data_to_test = get_available_meta_data(taxon_table_9_path, path_to_outdirs)

                    if meta_data_to_test == False:
                            sg.PopupError("Metadata table not found for the loaded TaXon table!\nPlease first create a metadata table!")
                    else:
                            win2_active = True
                            window.Hide()
                            meta_data_to_test_list = list(slices([sg.Radio(name, "metadata", key=name, default=True) for name in sorted(meta_data_to_test)], 5))

                            layout2 = [[sg.Text("NMDS analysis", size=(20,1))],
                                    [sg.Frame(layout = meta_data_to_test_list, title = 'Check taxa to exlude')],
                                    [sg.Button('Calculate')],
                                    [sg.Button('Back')]] # test = True

                            win2 = sg.Window('NMDS analysis', layout2, keep_on_top=False)

                            while True:
                                event2, values2 = win2.Read()
                                if event2 == 'Calculate':
                                        for input_value, test in values2.items():
                                                if test == True:
                                                        win2.Close()
                                                        win2_active = False
                                                        meta_data_to_test = input_value
                                                        from taxontabletools.NMDS_analysis import NMDS_analysis
                                                        NMDS_analysis(taxon_table_9_path, meta_data_to_test, nmds_analysis_taxonomic_level, width_nmds, height_nmds, nmds_s, max_iter_val, n_init_val, path_to_outdirs, template, font_size, color_discrete_sequence, nmds_dissimilarity)
                                                        window.UnHide()
                                                        break

                                if event2 is None or event2 == 'Back':
                                        win2.Close()
                                        win2_active = False
                                        window.UnHide()
                                        break

            if event == 'run_replicate_consistency_filter':
                if taxon_table_1_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    replicate_suffix_list = replicate_suffix_entry.replace(" ", "").split(",")
                    from taxontabletools.replicate_consistency_filter import replicate_consistency_filter
                    replicate_consistency_filter(taxon_table_1_path, replicate_suffix_list, path_to_outdirs, consistency_option)

            if event == 'run_create_perlodes_input_file':
                if taxon_table_11_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.convert_to_perlodes import convert_to_perlodes
                    convert_to_perlodes(taxon_table_11_path, operational_taxon_list_path, path_to_outdirs)

            if event == 'run_alpha_diversity':

                if taxon_table_8_path == '':
                    sg.PopupError("Please provide a file")
                else:

                    if alpha_diversity_taxonomic_level == "OTUs":
                        alpha_diversity_taxonomic_level = clustering_unit

                    from taxontabletools.get_available_meta_data import get_available_meta_data
                    meta_data_to_test = get_available_meta_data(taxon_table_8_path, path_to_outdirs)

                    if meta_data_to_test == False:
                            sg.PopupError("Metadata table not found for the loaded TaXon table!\nPlease first create a metadata table!")
                    else:
                            win2_active = True
                            window.Hide()
                            meta_data_to_test_list = list(slices([sg.Radio(name, "metadata", key=name, default=True) for name in sorted(meta_data_to_test)], 5))
                            layout2 = [[sg.Text("Alpha diversity", size=(20,1))],
                            [sg.Frame(layout = meta_data_to_test_list, title = 'Check metadata to test')],
                            [sg.Button('Calculate')],
                            [sg.Button('Back')]]

                            win2 = sg.Window('Alpha diversity estimates', layout2, keep_on_top=False)

                            while True:
                                event2, values2 = win2.Read()
                                if event2 == 'Calculate':
                                        for meta_data_to_test, test in values2.items():
                                            if test == True:
                                                if alpha_diversity_plot == "Box plot":
                                                        from taxontabletools.alpha_diversity import alpha_diversity_boxplot
                                                        alpha_diversity_boxplot(taxon_table_8_path, meta_data_to_test, width_alpha, height_alpha, alpha_diversity_taxonomic_level, path_to_outdirs, template, theme, font_size, color_discrete_sequence)
                                                        win2.Close()
                                                        win2_active = False
                                                        window.UnHide()
                                                        break
                                                if alpha_diversity_plot == "Scatter plot":
                                                        from taxontabletools.alpha_diversity import alpha_diversity_scatter_plot
                                                        alpha_diversity_scatter_plot(taxon_table_8_path, meta_data_to_test, width_alpha, height_alpha, alpha_scatter_size, alpha_diversity_taxonomic_level, path_to_outdirs, template, theme, font_size, color_discrete_sequence)
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
                    if taxon_table_8_path == '':
                        sg.PopupError("Please provide a file")
                    else:

                        if betadiv_heatmap_taxonomic_level == "OTUs":
                            betadiv_heatmap_taxonomic_level = clustering_unit

                        from taxontabletools.get_available_meta_data import get_available_meta_data
                        meta_data_to_test = get_available_meta_data(taxon_table_8_path, path_to_outdirs)

                        if meta_data_to_test == False:
                                sg.PopupError("Metadata table not found for the loaded TaXon table!\nPlease first create a metadata table!")
                        else:
                                win2_active = True
                                window.Hide()
                                meta_data_to_test_list = list(slices([sg.Radio(name, "metadata", key=name, default=True) for name in sorted(meta_data_to_test)], 5))

                                layout2 = [[sg.Text("Beta diversity", size=(20,1))],
                                [sg.Frame(layout = meta_data_to_test_list, title = 'Check metadata to test')],
                                [sg.Button('Calculate')],
                                [sg.Button('Back')]]

                                win2 = sg.Window('Beta diversity estimates', layout2, keep_on_top=False)

                                while True:
                                    event2, values2 = win2.Read()
                                    if event2 == 'Calculate':
                                            for meta_data_to_test, test in values2.items():
                                                if test == True:
                                                    from taxontabletools.beta_diversity import beta_diversity
                                                    beta_diversity(taxon_table_8_path, width_beta, height_beta, colorscale, meta_data_to_test, betadiv_heatmap_taxonomic_level, path_to_outdirs, template, font_size, betadiv_heatmap_dissimilarity)
                                                    win2.Close()
                                                    win2_active = False
                                                    window.UnHide()
                                                    break
                                    if event2 is None or event2 == 'Back':
                                            win2.Close()
                                            win2_active = False
                                            window.UnHide()
                                            break

            if event == 'run_betadiv_clustering':
                if taxon_table_8_path == '':
                    sg.PopupError("Please provide a file")
                else:
                    if betadiv_cluster_taxonomic_level == "OTUs":
                        betadiv_cluster_taxonomic_level = clustering_unit
                    from taxontabletools.beta_diversity import betadiv_clustering
                    betadiv_clustering(taxon_table_8_path, height_betadiv_clustering, width_betadiv_clustering, betadiv_clustering_threshold, betadiv_linkage, betadiv_cluster_taxonomic_level, path_to_outdirs, template, font_size, betadiv_cluster_dissimilarity)

            if event == 'run_parcat_analysis':
                if taxon_table_6_path == '':
                    sg.PopupError("Please provide a file")
                else:
                    if parcat_taxonomic_level == "OTUs":
                        parcat_taxonomic_level = clustering_unit

                    from taxontabletools.get_available_meta_data import get_available_meta_data
                    meta_data_to_test = get_available_meta_data(taxon_table_6_path, path_to_outdirs)

                    if meta_data_to_test == False:
                            sg.PopupError("Metadata table not found for the loaded TaXon table!\nPlease first create a metadata table!")
                    else:
                            win2_active = True
                            window.Hide()
                            meta_data_to_test_list = list(slices([sg.Radio(name, "metadata", key=name, default=True) for name in sorted(meta_data_to_test)], 5))
                            layout2 = [[sg.Text("Parallel category plot", size=(20,1))],
                            [sg.Frame(layout = meta_data_to_test_list, title = 'Check metadata to test')],
                            [sg.Button('Calculate')],
                            [sg.Button('Back')]]

                            win2 = sg.Window('Parallel category plot', layout2, keep_on_top=False)

                            while True:
                                event2, values2 = win2.Read()

                                if event2 == 'Calculate':
                                        for meta_data_to_test, test in values2.items():
                                            if test == True:
                                                win2.Close()
                                                win2_active = False
                                                from taxontabletools.parcat_analysis import parcat_analysis
                                                parcat_analysis(taxon_table_6_path, path_to_outdirs, template, height_parcat, width_parcat, meta_data_to_test, plotly_colors, available_taxonomic_levels_list, parcat_taxonomic_level, theme, font_size, color_discrete_sequence)
                                                window.UnHide()
                                if event2 is None or event2 == 'Back':
                                        win2.Close()
                                        win2_active = False
                                        window.UnHide()
                                        break

            if event == 'run_log_network':
                from taxontabletools.create_log import ttt_log_network
                ttt_log_network(path_to_outdirs)

            if event == 'run_open_log_file':
                try:
                    open_taxon_table(Path(str(path_to_outdirs) + "/log.xlsx"))
                except:
                    print(Path(str(path_to_outdirs) + "/log.xlsx"))
                    sg.PopupError("Log file not found!")

            if event == 'open_fgbewertung':
                webbrowser.open('https://www.gewaesser-bewertung-berechnung.de')

            if event == 'open_operational_taxonlist_download':
                webbrowser.open('https://www.gewaesser-bewertung-berechnung.de/index.php/perlodes-online.html')

            if event == 'open_github':
                webbrowser.open('https://github.com/TillMacher/TaxonTableTools')

            if event == 'open_twitter':
                webbrowser.open('https://twitter.com/THM_93')

            if event == 'open_tutorial':
                webbrowser.open('https://github.com/TillMacher/TaxonTableTools/blob/master/_data/TaxonTableTools_tutorial.pdf')

            if event == 'run_metadata_table_help':
                sg.Popup(metadata_table_help_text(), title = "Metadata table")

            if event == 'run_table_per_sample_help':
                sg.Popup(table_per_sample_help_text(), title = "Taxon Table per sample")

            if event == 'run_pa_conversion_help':
                sg.Popup(pa_conversion_help_text(), title = "Presence absence conversion")

            if event == 'run_alpha_diversity_help_text':
                sg.Popup(alpha_diversity_help_text(), title = "Alpha diversity")

            if event == 'run_nmds_analysis_help_text':
                sg.Popup(nmds_analysis_help_text(), title = "NMDS analysis")

            if event == 'run_beta_diversity_help_text':
                sg.Popup(beta_diversity_help_text(), title = "Beta diversity")

            if event == 'run_cluster_dendrogram_help_text':
                sg.Popup(cluster_dendrogram_help_text(), title = "Dendrogram clustering")

            if event == 'run_cca_analysis_help_text':
                sg.Popup(cca_analysis_help_text(), title = "CCA analysis")

            if event == 'run_pcoa_analysis_help_text':
                sg.Popup(pcoa_analysis_help_text(), title = "PCoA analysis")

            if event == 'run_layout_help':
                sg.Popup(layout_help_text(), title = "TaxonTableTools layout")

            if event == 'run_gbif_check_help':
                sg.Popup(gbif_check_help_text(), title = "GBIF check")

            if event == 'run_read_proportions_help':
                sg.Popup(read_proportions_help_text(), title = "How to optimize read proportion plots")

            if event == 'run_per_taxon_analysis':
                if taxon_table_6_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.per_taxon_analysis import per_taxon_analysis
                    per_taxon_analysis(taxon_table_6_path, height_per_taxon_stats, width_per_taxon_stats, per_taxon_stats_taxonomic_level, path_to_outdirs, template, theme, font_size, clustering_unit)

            if event == "run_read_filter":
                if taxon_table_2_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    if absolute_filtering == True:
                        read_filter_method = "absolute_filtering"
                        read_filter_treshold = absolute_filter_threshold
                    elif relative_filtering == True:
                        read_filter_method = "relative_filtering"
                        read_filter_treshold = relative_filter_threshold

                    from taxontabletools.taxon_table_filtering import read_filter
                    read_filter(taxon_table_2_path, path_to_outdirs, read_filter_method, read_filter_treshold)

            if event == "run_species_occurrence":
                if taxon_table_6_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    continents_to_check = []
                    if  africa_occurrence == True:
                        continents_to_check.append("Africa")
                    if  antarctica_occurrence == True:
                        continents_to_check.append("Antarctica")
                    if  asia_occurrence == True:
                        continents_to_check.append("Asia")
                    if  europe_occurrence == True:
                        continents_to_check.append("Europe")
                    if  north_america_occurrence == True:
                        continents_to_check.append("North America")
                    if  oceania_occurrence == True:
                        continents_to_check.append("Oceania")
                    if  south_america_occurrence == True:
                        continents_to_check.append("South America")

                    if continents_to_check != []:
                        from taxontabletools.gbif_occurrence import gbif_occurrence
                        gbif_occurrence(taxon_table_6_path, width_occurrence, height_occurrence, continents_to_check, template, theme, font_size, path_to_outdirs)
                    else:
                        sg.Popup("Please choose at least one continent", title="Error")

            if event == "run_species_occurrence_plot":
                from taxontabletools.gbif_occurrence import gbif_occurrence_plot
                gbif_occurrence_plot(width_occurrence, height_occurrence, template, theme, font_size, path_to_outdirs)

            if event == "run_show_colorsequence":
                fig = px.colors.qualitative.swatches()
                fig.show()

            if event == "run_show_colorscale":
                fig = px.colors.sequential.swatches()
                fig.show()

            if event == "open0":
                if taxon_table_0_path == '':
                    sg.PopupError("Please provide a file")
                else:
                    open_taxon_table(taxon_table_0_path)

            if event == "open1":
                if taxon_table_1_path == '':
                    sg.PopupError("Please provide a file")
                else:
                    open_taxon_table(taxon_table_1_path)

            if event == "open2":
                if taxon_table_2_path == '':
                    sg.PopupError("Please provide a file")
                else:
                    open_taxon_table(taxon_table_2_path)

            if event == "open3":
                if taxon_table_3_path == '':
                    sg.PopupError("Please provide a file")
                else:
                    open_taxon_table(taxon_table_3_path)

            if event == "open4":
                if taxon_table_4_path == '':
                    sg.PopupError("Please provide a file")
                else:
                    open_taxon_table(taxon_table_4_path)

            if event == "open5":
                if taxon_table_5_path == '':
                    sg.PopupError("Please provide a file")
                else:
                    open_taxon_table(taxon_table_5_path)

            if event == "open6":
                if taxon_table_6_path == '':
                    sg.PopupError("Please provide a file")
                else:
                    open_taxon_table(taxon_table_6_path)

            if event == "open7":
                if taxon_table_7_path == '':
                    sg.PopupError("Please provide a file")
                else:
                    open_taxon_table(taxon_table_7_path)

            if event == "open8":
                if taxon_table_8_path == '':
                    sg.PopupError("Please provide a file")
                else:
                    open_taxon_table(taxon_table_8_path)

            if event == "open9":
                if taxon_table_9_path == '':
                    sg.PopupError("Please provide a file")
                else:
                    open_taxon_table(taxon_table_9_path)

            if event == "open10":
                if taxon_table_10_path == '':
                    sg.PopupError("Please provide a file")
                else:
                    open_taxon_table(taxon_table_10_path)

            if event == "open11":
                if taxon_table_11_path == '':
                    sg.PopupError("Please provide a file")
                else:
                    open_taxon_table(taxon_table_11_path)

            if event == "run_save_preferences":
                f = open(user_preferences, "w")
                f.write(template + "\n")
                f.write(theme[0] + "\n")
                f.write(theme[1] + "\n")
                f.write(values["colorsequence"] + "\n")
                f.write(colorscale + "\n")
                f.write(str(font_size) + "\n")
                f.write(str(theme[2]) + "\n")
                f.write(clustering_unit)
                f.close()
                sg.Popup("Setting were saved...!")

            # close all windows to prevent crazy overlaps
            # only counts for matplotlib, but keep as long as venn is matplotlib-based
            if plt.get_fignums():
                plt.close()

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

        # Raised when an error is detected that doesn’t fall in any of the other categories. The associated value is a string indicating what precisely went wrong.
        # CHANGE TO CUSTOM ERROR and not RuntimeError!!!
        except RuntimeError:
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

        ###############################################################################################

    window.Close()

## run only if called as toplevel script
if __name__ == "__main__":
    main()
