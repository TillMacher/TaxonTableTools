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
taxon_tools_version = "Version 1.1.1"

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
    and is based on the number of OTUs per sample, which are displayed as scatter or bar plot.

    """
    return alpha_diversity_help_text

def beta_diversity_help_text():
    beta_diversity_help_text = """

    Beta diversity is calculated with the python package scikit-bio (http://scikit-bio.org/)
    and is calculated as jaccard-distances, which are illustrated in a distance matrix.
    An R-square test for reliability and validity is conducted.
    Check out matplotlib colormaps to customize the plot.

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
    Each sample from the selected TaXon table will be written to its own respective table.
    ____________________________________________________________________________

    This tool is particularly useful to compare single samples in the Venn diagram tool.
    OTUs that are not present in the respective sample are discarded.

    """
    return table_per_sample_help_text

def read_proportions_help_text():
    read_proportions_help_text = """
    Read proportion plots quickly become crowdy, unpleasing to view and hard to interprete.
    ____________________________________________________________________________

    The easiest workaround is to reduce the amount of data that is presented. Follow these
    few steps and the read proportions plots will quickly become much more eye-pleasing!
    1. Choosing a higher taxonomic taxonomic level
        -> On lower taxonomic levels (e.g. species) too many taxa are presented at once!
    2. Filter the data
        -> Removing of not required samples
        -> Removing of for the study irrelevant taxa groups
        -> Keep only species level hits (see Tutorial how to do that)

    Follow Ockham's razor principle as a guideline:
    -- The simplest explanation is usually the right one! --

    """
    return read_proportions_help_text

##########################################################################################################################
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

    directories_to_create = ["Venn_diagrams","TaXon_tables", "TaXon_tables_per_sample", "Taxon_lists", "Rarefaction_curves",
    "Taxonomic_richness_plots", "CCA_plots", "Taxonomic_resolution_plots", "Meta_data_table", "Site_occupancy_plots", "Read_proportions_plots", "Basic_stats",
    "Krona_charts", "Perlodes", "Alpha_diversity", "Beta_diversity", "PCoA_plots", "Replicate_analysis", "GBIF", "Occurrence_analysis", "Per_taxon_statistics"]

    for directory in directories_to_create:
        dirName = Path(str(path_to_outdirs) + "/" + directory + "/")
        if not os.path.exists(dirName):
            os.mkdir(Path(dirName))

    available_templates_list = ['seaborn', 'ggplot2', 'simple_white', 'plotly', 'plotly_dark', 'presentation', 'plotly_white']
    available_templates_radio = list(slices([sg.Radio(name, "template_radio", key=name, default=True) for name in available_templates_list], 7))
    available_color1_list = ["Dark", "Light", "Green", "Python", "Teal", "Blue"]
    available_color1_radio = list(slices([sg.Radio(name, "theme_radio", key=name, default=True) for name in available_color1_list], 6))

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

    ##########################################################################################################################
    ##########################################################################################################################
    ##########################################################################################################################

    getting_started_layout = [
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
    					[sg.Text('')],
    					[sg.Text('Graphic layout options', size=(50,2), font=('Arial', 11, "bold"))],
                        [sg.Frame(layout = available_templates_radio, title = 'Choose template')],
                        [sg.Frame(layout = available_color1_radio, title = 'Choose theme'), sg.CB("Transparent", key="transparent_theme", default=False)],
                        [sg.Frame(layout=[[sg.Text("Color 1:"), sg.Input("", key="custom_color1", size=(10,1)), sg.Text("Color 2:"), sg.Input("", key="custom_color2", size=(10,1)), sg.Radio("Custom", "theme_radio", key="custom_color")]], title="Or choose custom colors")],
    					[sg.Text('')],
    					]

    replicate_analysis_layout = [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Replicate processing', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'taxon_table_processing_1_table')],
    					[sg.Text('_'*105)],
    					[sg.Text('Replicates',size=(40,2), font=('Arial', 11, "bold"))],
                        [sg.Text("Replicate suffixes:", size=(30, 1)), sg.Input("a, b", size=(15,1), key="replicate_suffix_entry"),
                        sg.Text("", size=(1,1)), sg.Button("Preview", key="run_preview_suffix"), sg.Text("", size=(1,1)),
                        sg.Frame(layout=[[sg.Text("Comma separated list of suffixes.")]], title="Information")],
                        [sg.Text('',size=(1,1))],
    					[sg.Text('Replicate analysis',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Calculate shared OTUs", size=(30,1)), sg.Button("Run", key = 'run_replicate_analysis'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("1000", size=(4,1), key="width_repicates"),
                        sg.Input("600", size=(4,1), key="height_repicates"), sg.Text("", size=(1,1))]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text("Calculate OTU and read correlation", size=(30,1)), sg.Button("Run", key = 'run_replicate_correlation'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("1000", size=(4,1), key="width_repcorr"),
                        sg.Input("500", size=(4,1), key="height_repcorr"), sg.Text("", size=(1,1)), sg.Text("Start axes at zero"), sg.CB("x", key="x_zero_repcorr"), sg.CB("y", key="y_zero_repcorr")]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Replicate processing',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Replicate consistency", size=(30, 1)), sg.Button("Run", key='run_replicate_consistency_filter'), sg.Text("", size=(1,1)),
                        sg.Frame(layout=[[sg.Text("Keep OTUs present in all replicates.")]], title="Information")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text("Replicate merging", size=(30, 1)), sg.Button("Run", key='run_combine_replicates'), sg.Text("", size=(1,1)),
                        sg.Frame(layout=[[sg.Text("Merge replicates to one entry.")]], title="Information")],
    					[sg.Text('',size=(1,1))],
    					]

    taxon_table_filtering_layout = [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('TaXon table filtering', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'taxon_table_processing_2_table')],
    					[sg.Text('_'*105)],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Taxon-based filter',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Taxonomic level: "), sg.Radio('Phylum', "tt_filter", default=True, key='filter_phylum'), sg.Radio('Class', "tt_filter", key='filter_class'),
    					sg.Radio('Order', "tt_filter", key='filter_order'), sg.Radio('Family', "tt_filter", key='filter_family'),
    					sg.Radio('Genus', "tt_filter", key='filter_genus'), sg.Radio('Species', "tt_filter", key='filter_species')],
    					[sg.Text("", size=(20, 1))],
    					[sg.Text("Taxon-based filter:", size=(20, 1)), sg.Button("Run", key = 'run_taxon_filtering'), sg.Text("", size=(2,1)), sg.Frame(layout=[[sg.Radio('all', "cm_taxon_based_filter", key='cm_taxon_based_filter_all', default=True),
                        sg.Radio('none', "cm_taxon_based_filter", key='cm_taxon_based_filter_none')]], title="Check marks"), sg.Frame(layout=[[sg.Radio('Keep', "taxon_filter_method", key='taxon_filter_method_keep'),
                        sg.Radio('Exclude', "taxon_filter_method", key='taxon_filter_method_exclude', default=True)]], title="Method"), sg.Text("Mask:"), sg.Input(size=(5,1), key="taxon_filter_method_criterion")],
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
    					]

    data_conversion_layout = [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Data conversion', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'taxon_table_processing_3_table')],
    					[sg.Text('_'*105)],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('TaXon table per sample', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Create TaXon tables per-sample", size=(25,1)), sg.Button("Run", key = 'run_taxon_table_per_sample'), sg.Text("", size=(1,1)),
                        sg.Button("What is it good for?", button_color=('black', 'white'), key="run_table_per_sample_help" )],
    					[sg.Text('',size=(1,1))],
                		[sg.Text('Presence absence table', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Convert to p/a data", size=(25,1)), sg.Button("Run", key = 'run_convert_to_pa'), sg.Text("", size=(1,1)),
                        sg.Button("What is it good for?", button_color=('black', 'white'), key="run_pa_conversion_help" )],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Metadata table', size=(50,2), font=('Arial', 11, "bold")), sg.Text("", size=(2,1))],
    					[sg.Text("Create metadata table", size=(25,1)), sg.Button("Run", key = 'run_create_meta_data_table'), sg.Text("", size=(1,1)),
                        sg.Button("How does it work?", button_color=('black', 'white'), key="run_metadata_table_help" )],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('GBIF taxonomy check', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Check TaXon table:", size=(25, 1)), sg.Button("Run", key = 'run_gbif_check_taxonomy'), sg.Text("", size=(1,1)),
                        sg.Button("How does it work?", button_color=('black', 'white'), key="run_gbif_check_help" )],
    					[sg.Text('',size=(1,1))],
    					]

    basic_statistics_layout = [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Basic statistics', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'statistics_1_taxon_table_path')],
    					[sg.Text('_'*105)],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Basic statistics',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Calculate basic statistics", size=(25,1)), sg.Button("Run", key = 'run_basic_stats'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("800", size=(4,1), key="width_bstats"),
                        sg.Input("900", size=(4,1), key="height_bstats"), sg.Text("", size=(1,1))]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Taxonomic resolution',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Taxonomic resolution plot", size=(25,1)), sg.Button("Run", key = 'run_taxonomic_resolution'), sg.Text("", size=(1,1)),
                        sg.Frame(layout=[
                        [sg.Text("Plot size (w,h, ylim):"), sg.Input("600", size=(4,1), key="width_tax_res"), sg.Input("600", size=(4,1), key="height_tax_res")],
                        [sg.Text("Plot type:"), sg.Radio('a', "tres_plot_type", key='tres_type_a', default=True), sg.Radio('b', "tres_plot_type", key='tres_type_b')]
                        ], title="Settings", size=(1,2))],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Taxonomic richness',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Taxonomic richness plot", size=(25,1)), sg.Button("Run", key = 'run_taxonomic_richness'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h,ylim):"), sg.Input("600", size=(4,1), key="width_tax_rich"),
                        sg.Input("600", size=(4,1), key="height_tax_rich")]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Per taxon statistics',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Reads, OTUs and species", size=(25,1)), sg.Button("Run", key = 'run_per_taxon_analysis'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h,ylim):"), sg.Input("900", size=(4,1), key="width_per_taxon_stats"),
                        sg.Input("600", size=(4,1), key="height_per_taxon_stats")]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					]

    sample_comparison_layout =  [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Sample comparison', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'statistics_2_taxon_table_path')],
    					[sg.Text('_'*105)],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Venn diagrams',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text('Compare to table(s):', size=(16, 1)), sg.Input("", size=(10,1)), sg.FileBrowse(key = 'venn_taxon_table_2_path'), sg.Input("", size=(10,1)), sg.FileBrowse(key = 'venn_taxon_table_3_path'),
                        sg.Text("", size=(1,1)), sg.Text("Output file:"), sg.Input("venn_diagram", size=(15,1), key="venn_diagram_name")],
    					[sg.Text("", size=(1,1))],
    					[sg.Text("Compare TaXon tables", size=(20,1)), sg.Button("venn2", key = 'run_venn2_diagram'), sg.Button("venn3", key = 'run_venn3_diagram'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Text("1000", size=(4,1)),
                        sg.Text("1000", size=(4,1)), sg.Text("", size=(1,1))]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Site occupancy', size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Calculate site occupancy", size=(25,1)), sg.Button("Run", key = 'run_site_occupancy'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h): "), sg.Input("1000", size=(4,1), key="width_site_occ"),
                        sg.Input("1000", size=(4,1), key="height_site_occ"), sg.Text("", size=(1,1))]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Rarefaction curve',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Sample-based rarefaction", size=(25,1)), sg.Text("OTUs: ", size=(8,1)), sg.Button("Run", key = 'run_rarefaction_curve_OTUs'), sg.Text("", size=(1,1)), sg.Text("Species: "), sg.Button("Run", key = 'run_rarefaction_curve_species'), sg.Text("", size=(1,1)),
                        sg.Frame(layout=[[sg.Input('1000', size=(6,1), key='sample_repetitions'), sg.Text('repetitions')]], title="Settings")],
                        [sg.Text("Read-based rarefaction", size=(25,1)), sg.Text("Reads: ", size=(8,1)), sg.Button("Run", key = 'run_rarefaction_curve_reads'), sg.Text("", size=(1,1)),
                        sg.Frame(layout=[[sg.Input('10', size=(6,1), key='read_repetitions'), sg.Text('repetitions'), sg.Text("Plot size (w,h): "), sg.Input("1000", size=(4,1), key="width_rarefaction_reads"),
                        sg.Input("1000", size=(4,1), key="height_rarefaction_reads"), sg.Text("", size=(1,1))]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					]

    read_proportions_layout =  [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Read proportions', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'statistics_3_taxon_table_path')],
    					[sg.Text('_'*105)],
    					[sg.Text('',size=(1,1))],

    					[sg.Text('Read proportions',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Taxonomic level: "), sg.Radio('OTUs', "ra_filter", key='ra_OTUs'), sg.Radio('Phylum', "ra_filter", key='ra_phylum'), sg.Radio('Class', "ra_filter", key='ra_class'),
    					sg.Radio('Order', "ra_filter", key='ra_order'), sg.Radio('Family', "ra_filter", key='ra_family'),
    					sg.Radio('Genus', "ra_filter", key='ra_genus'), sg.Radio('Species', "ra_filter", default=True, key='ra_species')],
    					[sg.Text("", size=(20, 1))],

    					[sg.Text("Calculate heatmap", size=(25,1)), sg.Button("Run", key = 'run_read_proportions_heatmap'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h): "), sg.Input("1000", size=(4,1), key="width_read_props_heatmap"),
                        sg.Input("1000", size=(4,1), key="heigth_read_props_heatmap"), sg.Text("", size=(1,1))]], title="Settings")],
                        [sg.Text("Calculate bar chart", size=(25,1)), sg.Button("Run", key = 'run_read_proportions_bar'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h): "), sg.Input("1000", size=(4,1), key="width_read_props_bar"),
                        sg.Input("1000", size=(4,1), key="heigth_read_props_bar"), sg.Text("", size=(1,1))]], title="Settings"), sg.Button("How to optimize read proportion plots", key="run_read_proportions_help", button_color=("Black", "White"))],
                        [sg.Text("Calculate pie chart", size=(25,1)), sg.Button("Run", key = 'run_read_proportions_pie'), sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h): "), sg.Input("1000", size=(4,1), key="width_read_props_pie"),
                        sg.Input("1000", size=(4,1), key="heigth_read_props_pie"), sg.Text("", size=(1,1))]], title="Settings")],
    					[sg.Text('',size=(1,1))],

    					[sg.Text('Krona charts',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Create Krona chart", size=(25,1)), sg.Button("Run", key = 'run_create_krona_chart'), sg.Text("", size=(1,1)),
                        sg.Frame(layout=[[sg.Text("Charts are based on read proportions. Krona tools is required.")]], title="Information")],
    					[sg.Text('',size=(1,1))],
    					]

    diversity_analysis_layout =  [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Diversity and Ordination Analyses', size=(50,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'statistics_4_taxon_table_path')],
    					[sg.Text('_'*105)],
    					[sg.Text('',size=(1,1))],

    					[sg.Text('Diversity analyses',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Calculate alpha diversity", size=(23,1)), sg.Button("Run", key = 'run_alpha_diversity'),
                        sg.Button("Help", key = 'run_alpha_diversity_help_text', button_color=('black', 'white')),
                        sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("1000", size=(4,1), key="width_alpha"),
                        sg.Input("1000", size=(4,1), key="height_alpha")],
                        [sg.Radio("Scatter plot", "alpha_plot", key="alpha_scatter_plot"), sg.Text("Scatter size:"), sg.Input("40", size=(3,1), key="alpha_scatter_size")],
                        [sg.Radio("Box plot", "alpha_plot", key="alpha_box_plot", default=True)]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text("Calculate beta diversity", size=(23,1)), sg.Button("Run", key = 'run_beta_diversity'),
                        sg.Button("Help", key = 'run_beta_diversity_help_text', button_color=('black', 'white')),
                        sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("1000", size=(4,1), key="width_beta"),
                        sg.Input("1000", size=(4,1), key="height_beta")],
                        [sg.Text("Cmap:"), sg.Input("beta_cmap", size=(15,1), key="beta_cmap")]], title="Settings")],
    					[sg.Text('',size=(1,1))],

    					[sg.Text('Ordination analyses',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text("Calculate CCA", size=(23,1)), sg.Button("Run", key = 'run_CCA_analysis'),
                        sg.Button("Help", key = 'run_cca_analysis_help_text', button_color=('black', 'white')),
                        sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("1000", size=(4,1), key="width_cca"),
                        sg.Input("1000", size=(4,1), key="height_cca"), sg.Text("Scatter size:"), sg.Input("10", size=(3,1), key="cca_scatter_size"), sg.CB("Draw mesh", default=False, key="draw_mesh_cca")]], title="Settings")],
    					[sg.Text('',size=(1,1))],
    					[sg.Text("Calculate PCoA", size=(23,1)), sg.Button("Run", key = 'run_PCoA_analysis'),
                        sg.Button("Help", key = 'run_pcoa_analysis_help_text', button_color=('black', 'white')),
                        sg.Text("", size=(1,1)), sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("1000", size=(4,1), key="width_pcoa"),
                        sg.Input("1000", size=(4,1), key="height_pcoa"), sg.Text("Scatter size:"), sg.Input("10", size=(3,1), key="scatter_size_pcoa"),
                        sg.CB("Draw mesh", default=False, key="draw_mesh_pcoa")]], title="Settings")],
    					[sg.Text('',size=(1,1))],
        				]

    taxon_list_layout = [
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Taxon list',size=(40,2), font=('Arial', 11, "bold"))],
    					[sg.Text('TaXon table:', size=(20, 1)), sg.Input(), sg.FileBrowse(key = 'taxon_list_taxon_table_path')],
    					[sg.Text('_'*105)],
    					[sg.Text('',size=(1,1))],
    					[sg.Text('Convert to taxon list',size=(40,2), font=('Arial', 11, "bold"))],
                        [sg.Text("Enter name of output file:"), sg.Input(project_folder + "_taxon_list", size=(48,1), key='taxon_list_output_file_name'), sg.Text(".txt + .xlsx")],
                        [sg.Text('GBIF link:', size=(20,1)), sg.CB("", key="create_gbif_link", default=True), sg.Frame(layout=[[sg.Text("Internet connection is required.")]], title="Information")],
                        [sg.Text('Intraspecific distances:', size=(20,1)), sg.CB("", key="calc_dist", default=True), sg.Frame(layout=[[sg.Text("This option can drastically increase run times.")]], title="Information")],
                        [sg.Text('Occupancy per sample:', size=(20,1)), sg.CB("", key="calc_occupancy", default=True), sg.Frame(layout=[[sg.Text("This option can drastically increase run times.")]], title="Information")],
                        [sg.Text('', size=(1,1))],
                        [sg.Text("Create taxon list"), sg.Button("Run", key = 'run_create_taxon_list')],
                        [sg.Text('', size=(1,1))],
    					[sg.Text('Species occurrence analysis',size=(40,2), font=('Arial', 11, "bold"))],
                        [sg.Text('Collect species occurrences from GBIF'), sg.Button('Download', key="run_species_occurrence"), sg.Button('Plot', key="run_species_occurrence_plot"),
                        sg.Frame(layout=[[sg.Text("Plot size (w,h):"), sg.Input("1000", size=(4,1), key="width_occurrence"),
                        sg.Input("1000", size=(4,1), key="height_occurrence")]], title="Settings")],
                        [sg.Frame(layout=[[sg.CB("Africa", key="africa_occurrence"), sg.CB("Antarctica", key="antarctica_occurrence"),
                        sg.CB("Europe", key="europe_occurrence"), sg.CB("Asia", key="asia_occurrence"),
                        sg.CB("North America", key="north_america_occurrence"), sg.CB("Oceania", key="oceania_occurrence"),
                        sg.CB("South America", key="south_america_occurrence")]], title="Continents to check")],
    					[sg.Text('',size=(1,1))]
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
    					[sg.Button(key = 'open_fgbewertung', button_color=('white', 'white'), image_filename=fließgewaesserbewertungde)],
                        ]

    layout = [  [sg.Image(ttt_logo), sg.Text("", size=(9,1)), sg.Text('Project:', font=('Arial', 12, "bold")), sg.Text(project_folder, font=('Arial', 12, "bold"))],
    			[sg.Text('',size=(1,1))],
    			[sg.TabGroup([[
                sg.Tab('Getting started\n', getting_started_layout),
                sg.Tab('Replicate\nanalysis', replicate_analysis_layout),
                sg.Tab('TaXon table\nfiltering', taxon_table_filtering_layout),
                sg.Tab('Data\nconversion', data_conversion_layout),
                sg.Tab('Basic\nstatistics', basic_statistics_layout),
                sg.Tab('Sample\ncomparison', sample_comparison_layout),
                sg.Tab('Read\nproportions', read_proportions_layout),
                sg.Tab('Diversity\nanalyses', diversity_analysis_layout),
                sg.Tab('Taxon list\nconversion', taxon_list_layout),
                sg.Tab('Water Framework\nDirective', create_perlodes_input_layout)]])],
    			[sg.Text('',size=(1,1))],
    			[sg.Exit(button_color=('black', 'red')), sg.Text("", size=(35,1)),
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
            # collect plotly template
            template = [template for template in available_templates_list if values[template] == True][0]
            theme = [theme for theme in available_color1_list + ["custom_color"] if values[theme] == True][0]
            transparent_theme = values['transparent_theme']
            if transparent_theme != True:
                opacity_value = 1
            else:
                opacity_value = 0.5
            if theme == 'Dark':
                theme = ["Grey", "Black", opacity_value]
            if theme == 'Light':
                theme = ["White", "Grey", opacity_value]
            if theme == 'Green':
                theme = ["darkcyan", "Black", opacity_value]
            if theme == 'Teal':
                theme = ["teal", "Black", opacity_value]
            if theme == 'Python':
                theme = ["rgb(54,108,153)", "rgb(255,204,59)", opacity_value]
            if theme == 'Blue':
                theme = ["royalblue", "Black", opacity_value]
            if theme == "custom_color":
                color1 = values["custom_color1"]
                color2 = values["custom_color2"]
                if (color1.lower() in plotly_colors and color2.lower() in plotly_colors):
                    theme = [color1, color2, opacity_value]
                else:
                    sg.Popup("Warning: Colors not available. Choosing default theme.")
                    theme = ["royalblue", "Black", opacity_value]

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
            taxon_table_processing_3_table = values['taxon_table_processing_3_table']
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
            perlodes_taxon_table_path = values["perlodes_taxon_table_path"]
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
            calc_occupancy = values['calc_occupancy']
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
            taXon_table_to_check = values["taXon_table_to_check"]
            read_table_format_TTT = values["read_table_format_TTT"]
            read_table_format_jamp = values["read_table_format_jamp"]
            read_table_format_qiime = values["read_table_format_qiime"]
            alpha_box_plot = values["alpha_box_plot"]
            alpha_scatter_plot = values["alpha_scatter_plot"]
            height_repicates = values["height_repicates"]
            width_repicates = values["width_repicates"]
            width_repcorr = values["width_repcorr"]
            height_repcorr = values["height_repcorr"]
            height_per_taxon_stats = values["height_per_taxon_stats"]
            width_per_taxon_stats = values["width_per_taxon_stats"]
            beta_cmap = values["beta_cmap"]
            absolute_filtering = values["absolute_filtering"]
            absolute_filter_threshold = values["absolute_filter_threshold"]
            relative_filtering = values["relative_filtering"]
            relative_filter_threshold = values["relative_filter_threshold"]
            height_rarefaction_reads = values["height_rarefaction_reads"]
            width_rarefaction_reads = values["width_rarefaction_reads"]
            y_zero_repcorr = values["y_zero_repcorr"]
            x_zero_repcorr = values["x_zero_repcorr"]
            draw_mesh_pcoa = values["draw_mesh_pcoa"]
            draw_mesh_cca = values["draw_mesh_cca"]
            africa_occurrence = values["africa_occurrence"]
            antarctica_occurrence = values["antarctica_occurrence"]
            asia_occurrence = values["asia_occurrence"]
            europe_occurrence = values["europe_occurrence"]
            north_america_occurrence = values["north_america_occurrence"]
            oceania_occurrence = values["oceania_occurrence"]
            south_america_occurrence = values["south_america_occurrence"]
            height_occurrence = values["height_occurrence"]
            width_occurrence = values["width_occurrence"]

            print("########", "\n")

            if event == 'run_taxon_converter':

                if (taxonomy_table_path == '' or Read_table_path == ''):
                    sg.PopupError("Please provide a file", keep_on_top=True)
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
                elif read_table_format_jamp == True:
                    sg.Popup("Not supported yet!")
                elif read_table_format_qiime == True:
                    sg.Popup("Not supported yet!")

            if event == 'run_taXon_table_check':
                if (taXon_table_to_check == ''):
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.check_taXon_table_format import check_taXon_table_format
                    check_taXon_table_format(taXon_table_to_check)

            if event == "run_taxon_table_per_sample":
                if taxon_table_processing_3_table == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.create_taxon_table_per_sample import create_taxon_table_per_sample
                    create_taxon_table_per_sample(taxon_table_processing_3_table, path_to_outdirs)

            if event == "run_convert_to_pa":
                if taxon_table_processing_3_table == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.convert_to_presence_absence import convert_to_presence_absence
                    convert_to_presence_absence(taxon_table_processing_3_table, path_to_outdirs)

            if event == 'run_replicate_analysis':
                if taxon_table_processing_1_table == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    replicate_suffix_list = replicate_suffix_entry.replace(" ", "").split(",")
                    from taxontabletools.replicate_analysis import replicate_analysis
                    replicate_analysis(taxon_table_processing_1_table, height_repicates, width_repicates, replicate_suffix_list, path_to_outdirs, template, theme)

            if event == 'run_replicate_correlation':
                if taxon_table_processing_1_table == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    repcorr_options = [x_zero_repcorr, y_zero_repcorr]
                    replicate_suffix_list = replicate_suffix_entry.replace(" ", "").split(",")
                    from taxontabletools.replicate_correlation_analysis import replicate_correlation_analysis
                    replicate_correlation_analysis(taxon_table_processing_1_table, replicate_suffix_list, path_to_outdirs, width_repcorr, height_repcorr, repcorr_options, template, theme)

            if event == 'run_venn2_diagram':
                if os.path.exists(str(path_to_outdirs) + "/Venn_diagrams/" + venn_diagram_name):
                    sg.PopupError("Error: Venn directory already exists", keep_on_top=True)
                elif (statistics_2_taxon_table_path == '' or venn_taxon_table_2_path == ''):
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.venn_diagram import venn_diagram
                    venn_taxon_table_3_path = False
                    venn_diagram(statistics_2_taxon_table_path, venn_taxon_table_2_path, venn_taxon_table_3_path, venn_diagram_name, path_to_outdirs)

            if event == 'run_venn3_diagram':
                if os.path.exists(str(path_to_outdirs) + "/Venn_diagrams/" + venn_diagram_name):
                    sg.PopupError("Error: Venn directory already exists", keep_on_top=True)
                elif (statistics_2_taxon_table_path == '' or venn_taxon_table_2_path == ''):
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.venn_diagram import venn_diagram
                    venn_diagram(statistics_2_taxon_table_path, venn_taxon_table_2_path, venn_taxon_table_3_path, venn_diagram_name, path_to_outdirs)

            if event == 'run_rarefaction_curve_OTUs':
                if statistics_2_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.rarefaction_curve import rarefaction_curve_OTUs
                    rarefaction_curve_OTUs(statistics_2_taxon_table_path, int(sample_repetitions), path_to_outdirs, template, theme)

            if event == 'run_rarefaction_curve_species':
                if statistics_2_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.rarefaction_curve import rarefaction_curve_species
                    rarefaction_curve_species(statistics_2_taxon_table_path, int(sample_repetitions), path_to_outdirs, template, theme)

            if event == 'run_rarefaction_curve_reads':
                if statistics_2_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    answer = sg.PopupOKCancel("This analysis can take very long. Continue?", title="Runtime")
                    if answer == "OK":
                        from taxontabletools.rarefaction_curve import rarefaction_curve_reads
                        rarefaction_curve_reads(statistics_2_taxon_table_path, int(read_repetitions), width_rarefaction_reads, height_rarefaction_reads, path_to_outdirs, template, theme)

            if event == 'run_create_taxon_list':
                if taxon_list_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.create_taxon_list import create_taxon_list
                    create_taxon_list(taxon_list_taxon_table_path, taxon_list_output_file_name, create_gbif_link, calc_dist, calc_occupancy, taxon_tools_version, path_to_outdirs)

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
                else:
                    from taxontabletools.taxon_table_filtering import create_taxa_mask
                    available_taxa = create_taxa_mask(taxon_table_processing_2_table, mask, taxon_filter_method_criterion)
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
                                taxon_filter(taxon_table_processing_2_table, check_list, mask, appendix_name, threshold, path_to_outdirs, taxon_filter_method)
                            else:
                                taxon_filter(taxon_table_processing_2_table, check_list, mask, appendix_name, threshold, path_to_outdirs, taxon_filter_method)

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
                else:
                    from taxontabletools.taxon_table_filtering import create_sample_mask
                    available_samples = sorted(create_sample_mask(taxon_table_processing_2_table, sample_filter_method_criterion))
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
                            filter_samples(taxon_table_processing_2_table, check_list, appendix_name, path_to_outdirs, sample_filter_method)

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
                if taxon_table_processing_3_table == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.gbif_check_taxonomy import gbif_check_taxonomy
                    gbif_check_taxonomy(taxon_table_processing_3_table, path_to_outdirs)

            if event == 'run_combine_replicates':

                replicate_suffix_list = replicate_suffix_entry.replace(" ", "").split(",")

                if taxon_table_processing_1_table == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
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
                else:
                    from taxontabletools.calculate_taxonomic_richness import calculate_taxonomic_richness
                    calculate_taxonomic_richness(statistics_1_taxon_table_path, path_to_outdirs, width_tax_rich, height_tax_rich, template, theme)

            if event == 'run_taxonomic_resolution':
                if statistics_1_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    if tres_type_a == True:
                        figure_type = "a"
                    if tres_type_b == True:
                        figure_type = "b"
                    from taxontabletools.calculate_taxonomic_resolution import calculate_taxonomic_resolution
                    calculate_taxonomic_resolution(statistics_1_taxon_table_path, path_to_outdirs, width_tax_res, height_tax_res, figure_type, template, theme)

            if event == 'run_create_meta_data_table':
                if taxon_table_processing_3_table == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.create_metadata_table import create_metadata_table
                    create_metadata_table(taxon_table_processing_3_table, path_to_outdirs)

            if event == 'run_site_occupancy':
                if statistics_2_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.get_available_meta_data import get_available_meta_data
                    meta_data_to_test = get_available_meta_data(statistics_2_taxon_table_path, path_to_outdirs)

                    if meta_data_to_test == False:
                                 sg.PopupError("Metadata table not found for the loaded TaXon table!\nPlease first create a metadata table!", keep_on_top=True)
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
                                                           site_occupancy(statistics_2_taxon_table_path, meta_data_to_test, taxon_level, path_to_outdirs, width_site_occ, height_site_occ, template, theme)
                                                           win2.Close()
                                                           win2_active = False
                                                           window.UnHide()

                                              if event2 is None or event2 == 'Back':
                                                           win2.Close()
                                                           win2_active = False
                                                           window.UnHide()
                                                           break

                                 print("\n" + "########")

            if event == 'run_read_proportions_heatmap':
                if statistics_3_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
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

                    from taxontabletools.read_proportions import read_proportions_heatmap
                    read_proportions_heatmap(statistics_3_taxon_table_path, taxonomic_level, path_to_outdirs, width_read_props_heatmap, heigth_read_props_heatmap, template)

            if event == 'run_read_proportions_bar':
                if statistics_3_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
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

                    from taxontabletools.read_proportions import read_proportions_bar
                    read_proportions_bar(statistics_3_taxon_table_path, taxonomic_level, path_to_outdirs, width_read_props_bar, heigth_read_props_bar, template)

            if event == 'run_read_proportions_pie':
                if statistics_3_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
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

                    from taxontabletools.read_proportions import read_proportions_pie
                    read_proportions_pie(statistics_3_taxon_table_path, taxonomic_level, path_to_outdirs, width_read_props_pie, heigth_read_props_pie, template)

            if event == 'run_basic_stats':
                if statistics_1_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.basic_stats import basic_stats
                    basic_stats(statistics_1_taxon_table_path, height_bstats, width_bstats, path_to_outdirs, template, theme)

            if event == 'run_create_krona_chart':
                if statistics_3_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.create_krona_chart import create_krona_chart
                    create_krona_chart(statistics_3_taxon_table_path, path_to_outdirs)

            if event == 'run_CCA_analysis':
                if statistics_4_taxon_table_path == '':
                    sg.PopupError("Please provide a file")
                else:
                    from taxontabletools.get_available_meta_data import get_available_meta_data
                    meta_data_to_test = get_available_meta_data(statistics_4_taxon_table_path, path_to_outdirs)

                    if meta_data_to_test == False:
                            sg.PopupError("Metadata table not found for the loaded TaXon table!\nPlease first create a metadata table!")
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
                                                        CCA_analysis(statistics_4_taxon_table_path, meta_data_to_test, width_cca, height_cca, cca_scatter_size, draw_mesh_cca, path_to_outdirs, template)
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
                else:
                    from taxontabletools.get_available_meta_data import get_available_meta_data
                    meta_data_to_test = get_available_meta_data(statistics_4_taxon_table_path, path_to_outdirs)

                    if meta_data_to_test == False:
                            sg.PopupError("Metadata table not found for the loaded TaXon table!\nPlease first create a metadata table!")
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
                                                        PCoA_analysis(statistics_4_taxon_table_path, meta_data_to_test, width_pcoa, height_pcoa, scatter_size_pcoa, draw_mesh_pcoa, path_to_outdirs, template)
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
                else:
                    replicate_suffix_list = replicate_suffix_entry.replace(" ", "").split(",")
                    from taxontabletools.replicate_consistency_filter import replicate_consistency_filter
                    replicate_consistency_filter(taxon_table_processing_1_table, replicate_suffix_list, path_to_outdirs)

            if event == 'run_create_perlodes_input_file':
                if perlodes_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    from taxontabletools.convert_to_perlodes import convert_to_perlodes
                    convert_to_perlodes(perlodes_taxon_table_path, operational_taxon_list_path, path_to_outdirs)

            if event == 'run_alpha_diversity':

                    if statistics_4_taxon_table_path == '':
                        sg.PopupError("Please provide a file")
                    else:
                        from taxontabletools.get_available_meta_data import get_available_meta_data
                        meta_data_to_test = get_available_meta_data(statistics_4_taxon_table_path, path_to_outdirs)

                        if meta_data_to_test == False:
                                sg.PopupError("Metadata table not found for the loaded TaXon table!\nPlease first create a metadata table!")
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
                                                            alpha_diversity_boxplot(statistics_4_taxon_table_path, meta_data_to_test, width_alpha, height_alpha, alpha_scatter_size, path_to_outdirs, template, theme)
                                                            win2.Close()
                                                            win2_active = False
                                                            window.UnHide()
                                                            break
                                                    if alpha_scatter_plot == True:
                                                            from taxontabletools.alpha_diversity import alpha_diversity_scatter_plot
                                                            alpha_diversity_scatter_plot(statistics_4_taxon_table_path, meta_data_to_test, width_alpha, height_alpha, alpha_scatter_size, path_to_outdirs, template)
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
                    else:
                        from taxontabletools.get_available_meta_data import get_available_meta_data
                        meta_data_to_test = get_available_meta_data(statistics_4_taxon_table_path, path_to_outdirs)

                        if meta_data_to_test == False:
                                sg.PopupError("Metadata table not found for the loaded TaXon table!\nPlease first create a metadata table!")
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
                                                    beta_diversity(statistics_4_taxon_table_path, width_beta, height_beta, beta_cmap, meta_data_to_test, path_to_outdirs, template)
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

            if event == 'run_beta_diversity_help_text':
                sg.Popup(beta_diversity_help_text(), title = "Beta diversity")

            if event == 'run_cca_analysis_help_text':
                sg.Popup(cca_analysis_help_text(), title = "CCA analysis")

            if event == 'run_pcoa_analysis_help_text':
                sg.Popup(pcoa_analysis_help_text(), title = "PCoA analysis")

            if event == 'run_gbif_check_help':
                sg.Popup(gbif_check_help_text(), title = "GBIF check")

            if event == 'run_read_proportions_help':
                sg.Popup(read_proportions_help_text(), title = "How to optimize read proportion plots")

            if event == 'run_per_taxon_analysis':
                if statistics_1_taxon_table_path == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    win2_active = True
                    window.Hide()
                    taxonmic_levels_to_choose = ['Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
                    taxonmic_levels_list = list(slices([sg.Radio(name, "taxa", key=name, default=True) for name in taxonmic_levels_to_choose], 6))
                    layout2 = [[sg.Text("Per taxon statistics", size=(20,1))],
                    [sg.Text('',size=(1,1))],
                    [sg.Frame(layout = taxonmic_levels_list, title = 'Check taxonomic level to use')],
                    [sg.Text('',size=(1,1))],
                    [sg.Button('Calculate')],
                    [sg.Button('Back')]] # test = True

                    win2 = sg.Window('Per taxon statistics', layout2, keep_on_top=True)

                    while True:
                        event2, values2 = win2.Read()
                        if event2 == 'Calculate':
                            for input_value, test in values2.items():
                                if (test == True and input_value in taxonmic_levels_to_choose):
                                    taxonomic_level = input_value
                            win2_active = False
                            window.UnHide()
                            win2.Close()
                            from taxontabletools.per_taxon_analysis import per_taxon_analysis
                            per_taxon_analysis(statistics_1_taxon_table_path, height_per_taxon_stats, width_per_taxon_stats, taxonomic_level, path_to_outdirs, template, theme)
                        if event2 is None or event2 == 'Back':
                            win2.Close()
                            win2_active = False
                            window.UnHide()
                            break

            if event == "run_read_filter":
                if taxon_table_processing_2_table == '':
                    sg.PopupError("Please provide a file", keep_on_top=True)
                else:
                    if absolute_filtering == True:
                        read_filter_method = "absolute_filtering"
                        read_filter_treshold = absolute_filter_threshold
                    elif relative_filtering == True:
                        read_filter_method = "relative_filtering"
                        read_filter_treshold = relative_filter_threshold

                    from taxontabletools.taxon_table_filtering import read_filter
                    read_filter(taxon_table_processing_2_table, path_to_outdirs, read_filter_method, read_filter_treshold)

            if event == "run_species_occurrence":
                if taxon_list_taxon_table_path == '':
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
                        gbif_occurrence(taxon_list_taxon_table_path, width_occurrence, height_occurrence, continents_to_check, template, theme, path_to_outdirs)
                    else:
                        sg.Popup("Please choose at least one continent", title="Error")

            if event == "run_species_occurrence_plot":
                from taxontabletools.gbif_occurrence import gbif_occurrence_plot
                gbif_occurrence_plot(width_occurrence, height_occurrence, template, theme, path_to_outdirs)

            # close all windows to prevent crazy overlaps
            # only counts for matplotlib, but keep as long as venn is matplotlib-based
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
