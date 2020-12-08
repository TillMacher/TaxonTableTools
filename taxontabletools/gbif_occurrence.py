def gbif_occurrence(TaXon_table_xlsx, width, height, continents_to_check, template, theme, path_to_outdirs):

    import requests_html, json
    import PySimpleGUI as sg
    import pandas as pd
    from pandas import DataFrame
    import numpy as np
    from pathlib import Path
    import plotly.graph_objects as go
    import os

    ## dictionary with all country codes of the Earth
    country_codes_dict = {'Andorra': ['AD', 'Europe'],
 'United Arab Emirates': ['AE', 'Asia'],
 'Afghanistan': ['AF', 'Asia'],
 'Antigua and Barbuda': ['AG', 'North America'],
 'Anguilla': ['AI', 'North America'],
 'Albania': ['AL', 'Europe'],
 'Armenia': ['AM', 'Asia'],
 'Angola': ['AO', 'Africa'],
 'Antarctica': ['AQ', 'Antarctica'],
 'Argentina': ['AR', 'South America'],
 'American Samoa': ['AS', 'Oceania'],
 'Austria': ['AT', 'Europe'],
 'Australia': ['AU', 'Oceania'],
 'Aruba': ['AW', 'North America'],
 'Åland Islands': ['AX', 'Europe'],
 'Azerbaijan': ['AZ', 'Asia'],
 'Bosnia and Herzegovina': ['BA', 'Europe'],
 'Barbados': ['BB', 'North America'],
 'Bangladesh': ['BD', 'Asia'],
 'Belgium': ['BE', 'Europe'],
 'Burkina Faso': ['BF', 'Africa'],
 'Bulgaria': ['BG', 'Europe'],
 'Bahrain': ['BH', 'Asia'],
 'Burundi': ['BI', 'Africa'],
 'Benin': ['BJ', 'Africa'],
 'Saint Barthélemy': ['BL', 'North America'],
 'Bermuda': ['BM', 'North America'],
 'Brunei Darussalam': ['BN', 'Asia'],
 'Bolivia': ['BO', 'South America'],
 'Bonaire, Sint Eustatius and Saba': ['BQ', 'North America'],
 'Brazil': ['BR', 'South America'],
 'Bahamas': ['BS', 'North America'],
 'Bhutan': ['BT', 'Asia'],
 'Bouvet Island': ['BV', 'Antarctica'],
 'Botswana': ['BW', 'Africa'],
 'Belarus': ['BY', 'Europe'],
 'Belize': ['BZ', 'North America'],
 'Canada': ['CA', 'North America'],
 'Cocos (Keeling) Islands': ['CC', 'Asia'],
 'Congo (Democratic Republic)': ['CD', 'Africa'],
 'Central African Republic': ['CF', 'Africa'],
 'Congo': ['CG', 'Africa'],
 'Switzerland': ['CH', 'Europe'],
 "Côte d'Ivoire": ['CI', 'Africa'],
 'Cook Islands': ['CK', 'Oceania'],
 'Chile': ['CL', 'South America'],
 'Cameroon': ['CM', 'Africa'],
 'China': ['CN', 'Asia'],
 'Colombia': ['CO', 'South America'],
 'Costa Rica': ['CR', 'North America'],
 'Cuba': ['CU', 'North America'],
 'Cabo Verde': ['CV', 'Africa'],
 'Curaçao': ['CW', 'North America'],
 'Christmas Island': ['CX', 'Asia'],
 'Cyprus': ['CY', 'Asia'],
 'Czechia': ['CZ', 'Europe'],
 'Germany': ['DE', 'Europe'],
 'Djibouti': ['DJ', 'Africa'],
 'Denmark': ['DK', 'Europe'],
 'Dominica': ['DM', 'North America'],
 'Dominican Republic': ['DO', 'North America'],
 'Algeria': ['DZ', 'Africa'],
 'Ecuador': ['EC', 'South America'],
 'Estonia': ['EE', 'Europe'],
 'Egypt': ['EG', 'Africa'],
 'Western Sahara': ['EH', 'Africa'],
 'Eritrea': ['ER', 'Africa'],
 'Spain': ['ES', 'Europe'],
 'Ethiopia': ['ET', 'Africa'],
 'Finland': ['FI', 'Europe'],
 'Fiji': ['FJ', 'Oceania'],
 'Falkland Islands': ['FK', 'South America'],
 'Micronesia': ['FM', 'Oceania'],
 'Faroe Islands': ['FO', 'Europe'],
 'France': ['FR', 'Europe'],
 'Gabon': ['GA', 'Africa'],
 'United Kingdom': ['GB',
  'Europe'],
 'Grenada': ['GD', 'North America'],
 'Georgia': ['GE', 'Asia'],
 'French Guiana': ['GF', 'South America'],
 'Guernsey': ['GG', 'Europe'],
 'Ghana': ['GH', 'Africa'],
 'Gibraltar': ['GI', 'Europe'],
 'Greenland': ['GL', 'North America'],
 'Gambia': ['GM', 'Africa'],
 'Guinea': ['GN', 'Africa'],
 'Guadeloupe': ['GP', 'North America'],
 'Equatorial Guinea': ['GQ', 'Africa'],
 'Greece': ['GR', 'Europe'],
 'South Georgia and the South Sandwich Islands': ['GS', 'Antarctica'],
 'Guatemala': ['GT', 'North America'],
 'Guam': ['GU', 'Oceania'],
 'Guinea-Bissau': ['GW', 'Africa'],
 'Guyana': ['GY', 'South America'],
 'Hong Kong': ['HK', 'Asia'],
 'Heard Island and McDonald Islands': ['HM', 'Antarctica'],
 'Honduras': ['HN', 'North America'],
 'Croatia': ['HR', 'Europe'],
 'Haiti': ['HT', 'North America'],
 'Hungary': ['HU', 'Europe'],
 'Indonesia': ['ID', 'Asia'],
 'Ireland': ['IE', 'Europe'],
 'Israel': ['IL', 'Asia'],
 'Isle of Man': ['IM', 'Europe'],
 'India': ['IN', 'Asia'],
 'British Indian Ocean Territory': ['IO', 'Asia'],
 'Iraq': ['IQ', 'Asia'],
 'Iran': ['IR', 'Asia'],
 'Iceland': ['IS', 'Europe'],
 'Italy': ['IT', 'Europe'],
 'Jersey': ['JE', 'Europe'],
 'Jamaica': ['JM', 'North America'],
 'Jordan': ['JO', 'Asia'],
 'Japan': ['JP', 'Asia'],
 'Kenya': ['KE', 'Africa'],
 'Kyrgyzstan': ['KG', 'Asia'],
 'Cambodia': ['KH', 'Asia'],
 'Kiribati': ['KI', 'Oceania'],
 'Comoros': ['KM', 'Africa'],
 'Saint Kitts and Nevis': ['KN', 'North America'],
 "Korea (Democratic People's Republic)": ['KP', 'Asia'],
 'Korea (Republic)': ['KR', 'Asia'],
 'Kuwait': ['KW', 'Asia'],
 'Cayman Islands': ['KY', 'North America'],
 'Kazakhstan': ['KZ', 'Asia'],
 "Lao People's Democratic Republic": ['LA', 'Asia'],
 'Lebanon': ['LB', 'Asia'],
 'Saint Lucia': ['LC', 'North America'],
 'Liechtenstein': ['LI', 'Europe'],
 'Sri Lanka': ['LK', 'Asia'],
 'Liberia': ['LR', 'Africa'],
 'Lesotho': ['LS', 'Africa'],
 'Lithuania': ['LT', 'Europe'],
 'Luxembourg': ['LU', 'Europe'],
 'Latvia': ['LV', 'Europe'],
 'Libya': ['LY', 'Africa'],
 'Morocco': ['MA', 'Africa'],
 'Monaco': ['MC', 'Europe'],
 'Moldova (the Republic of)': ['MD', 'Europe'],
 'Montenegro': ['ME', 'Europe'],
 'Saint Martin (French part)': ['MF', 'North America'],
 'Madagascar': ['MG', 'Africa'],
 'Marshall Islands': ['MH', 'Oceania'],
 'Republic of North Macedonia': ['MK', 'Europe'],
 'Mali': ['ML', 'Africa'],
 'Myanmar': ['MM', 'Asia'],
 'Mongolia': ['MN', 'Asia'],
 'Macao': ['MO', 'Asia'],
 'Northern Mariana Islands': ['MP', 'Oceania'],
 'Martinique': ['MQ', 'North America'],
 'Mauritania': ['MR', 'Africa'],
 'Montserrat': ['MS', 'North America'],
 'Malta': ['MT', 'Europe'],
 'Mauritius': ['MU', 'Africa'],
 'Maldives': ['MV', 'Asia'],
 'Malawi': ['MW', 'Africa'],
 'Mexico': ['MX', 'North America'],
 'Malaysia': ['MY', 'Asia'],
 'Mozambique': ['MZ', 'Africa'],
 'Namibia': ['NA', 'Africa'],
 'New Caledonia': ['NC', 'Oceania'],
 'Niger': ['NE', 'Africa'],
 'Norfolk Island': ['NF', 'Oceania'],
 'Nigeria': ['NG', 'Africa'],
 'Nicaragua': ['NI', 'North America'],
 'Netherlands': ['NL', 'Europe'],
 'Norway': ['NO', 'Europe'],
 'Nepal': ['NP', 'Asia'],
 'Nauru': ['NR', 'Oceania'],
 'Niue': ['NU', 'Oceania'],
 'New Zealand': ['NZ', 'Oceania'],
 'Oman': ['OM', 'Asia'],
 'Panama': ['PA', 'North America'],
 'Peru': ['PE', 'South America'],
 'French Polynesia': ['PF', 'Oceania'],
 'Papua New Guinea': ['PG', 'Oceania'],
 'Philippines': ['PH', 'Asia'],
 'Pakistan': ['PK', 'Asia'],
 'Poland': ['PL', 'Europe'],
 'Saint Pierre and Miquelon': ['PM', 'North America'],
 'Pitcairn': ['PN', 'Oceania'],
 'Puerto Rico': ['PR', 'North America'],
 'Palestine, State of': ['PS', 'Asia'],
 'Portugal': ['PT', 'Europe'],
 'Palau': ['PW', 'Oceania'],
 'Paraguay': ['PY', 'South America'],
 'Qatar': ['QA', 'Asia'],
 'Réunion': ['RE', 'Africa'],
 'Romania': ['RO', 'Europe'],
 'Serbia': ['RS', 'Europe'],
 'Russian Federation': ['RU', 'Europe'],
 'Rwanda': ['RW', 'Africa'],
 'Saudi Arabia': ['SA', 'Asia'],
 'Solomon Islands': ['SB', 'Oceania'],
 'Seychelles': ['SC', 'Africa'],
 'Sudan': ['SD', 'Africa'],
 'Sweden': ['SE', 'Europe'],
 'Singapore': ['SG', 'Asia'],
 'Saint Helena, Ascension and Tristan da Cunha': ['SH', 'Africa'],
 'Slovenia': ['SI', 'Europe'],
 'Svalbard and Jan Mayen': ['SJ', 'Europe'],
 'Slovakia': ['SK', 'Europe'],
 'Sierra Leone': ['SL', 'Africa'],
 'San Marino': ['SM', 'Europe'],
 'Senegal': ['SN', 'Africa'],
 'Somalia': ['SO', 'Africa'],
 'Suriname': ['SR', 'South America'],
 'South Sudan': ['SS', 'Africa'],
 'Sao Tome and Principe': ['ST', 'Africa'],
 'El Salvador': ['SV', 'North America'],
 'Syrian Arab Republic': ['SY', 'Asia'],
 'Eswatini': ['SZ', 'Africa'],
 'Turks and Caicos Islands': ['TC', 'North America'],
 'Chad': ['TD', 'Africa'],
 'French Southern Territories': ['TF', 'Antarctica'],
 'Togo': ['TG', 'Africa'],
 'Thailand': ['TH', 'Asia'],
 'Tajikistan': ['TJ', 'Asia'],
 'Tokelau': ['TK', 'Oceania'],
 'Timor-Leste': ['TL', 'Asia'],
 'Turkmenistan': ['TM', 'Asia'],
 'Tunisia': ['TN', 'Africa'],
 'Tonga': ['TO', 'Oceania'],
 'Turkey': ['TR', 'Europe'],
 'Trinidad and Tobago': ['TT', 'North America'],
 'Tuvalu': ['TV', 'Oceania'],
 'Taiwan': ['TW', 'Asia'],
 'Tanzania': ['TZ', 'Africa'],
 'Ukraine': ['UA', 'Europe'],
 'Uganda': ['UG', 'Africa'],
 'United States Minor Outlying Islands': ['UM', 'Oceania'],
 'United States of America': ['US', 'North America'],
 'Uruguay': ['UY', 'South America'],
 'Uzbekistan': ['UZ', 'Asia'],
 'Holy See': ['VA', 'Europe'],
 'Saint Vincent and the Grenadines': ['VC', 'North America'],
 'Venezuela (Bolivarian Republic of)': ['VE', 'South America'],
 'Virgin Islands (British)': ['VG', 'North America'],
 'Virgin Islands (U.S.)': ['VI', 'North America'],
 'Viet Nam': ['VN', 'Asia'],
 'Vanuatu': ['VU', 'Oceania'],
 'Wallis and Futuna': ['WF', 'Oceania'],
 'Samoa': ['WS', 'Oceania'],
 'Yemen': ['YE', 'Asia'],
 'Mayotte': ['YT', 'Africa'],
 'South Africa': ['ZA', 'Africa'],
 'Zambia': ['ZM', 'Africa'],
 'Zimbabwe': ['ZW', 'Africa']}

    ## load Taxon table
    #TaXon_table_xlsx = "/Users/tillmacher/Desktop/Projects/TTT_Projects/Projects/Sicliy_MZB/TaXon_tables/Sicily_eDNA_MZB_taxon_table_renamed_cons_derep_no_match_excluded_blanks_excluded_species.xlsx"
    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx).fillna("nan")
    ## get a unique list of taxa
    taxa = set([taxon for taxon in TaXon_table_df["Species"].values.tolist() if taxon != "nan"])
    ## create a dataframe to store the results
    ## select only countries for the selected continents
    selected_countries_list = [country for country, values in country_codes_dict.items() if values[1] in continents_to_check]
    occurrence_df = pd.DataFrame(selected_countries_list, columns=["Country"])

    ## calculate runtime
    n_countries = len(selected_countries_list)
    n_species = len(taxa)
    t_single_request = 0.15
    t_total = round(t_single_request * n_species * n_countries / 60, 1)
    ## ask to continue the script
    answer = sg.PopupOKCancel("This will take roughly " + str(t_total) + " minutes. Continue?", title="Runtime")

    if answer == 'OK':

        ## create a subfolder for better sorting and overview
        dirName = Path(str(path_to_outdirs) + "/" + "Occurrence_analysis" + "/" + TaXon_table_xlsx.stem + "/")
        if not os.path.exists(dirName):
            os.mkdir(dirName)

        ############################################################################
        ## create the progress bar window
        layout = [[sg.Text('Progress bar')],
                  [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
                  [sg.Cancel()]]
        window_progress_bar = sg.Window('Progress bar', layout, keep_on_top=True)
        progress_bar = window_progress_bar['progressbar']
        progress_update = 0
        progress_increase = 1000 / len(taxa) + 1
        ############################################################################

        ############################################################################
        event, values = window_progress_bar.read(timeout=10)
        if event == 'Cancel'  or event is None:
            window_progress_bar.Close()
            raise RuntimeError
        # update bar with loop value +1 so that bar eventually reaches the maximum
        progress_update += 0
        progress_bar.UpdateBar(progress_update)
        ############################################################################

        ## request gbif for the occurrence data
        for taxon_name in taxa:
            occurrence_list= []
            for country, values in country_codes_dict.items():
                country_code = values[0]
                continent = values[1]
                ## only check selected continents to reduce runtime
                if continent in continents_to_check:
                    ## create an html session
                    with requests_html.HTMLSession() as session:
                        ## generate html request name
                        request_name = '%20'.join(taxon_name.split(' '))
                        ## request that name
                        r = session.get("https://api.gbif.org/v1/occurrence/search?scientificName=" + request_name + "&country=" + country_code)
                        ## parse json
                        res = json.loads(r.text)
                        ## get number of occurrences
                        occurrence_list.append(res["count"])
            ## store the results in the dataframe
            occurrence_df[taxon_name] = occurrence_list

            ############################################################################
            event, values = window_progress_bar.read(timeout=10)
            if event == 'Cancel'  or event is None:
                window_progress_bar.Close()
                raise RuntimeError
            # update bar with loop value +1 so that bar eventually reaches the maximum
            progress_update += progress_increase
            progress_bar.UpdateBar(progress_update)
            ############################################################################

        window_progress_bar.Close()

        ## remove countries that have 0 hits
        occurrence_df_filtered_list = []
        for row in occurrence_df.values.tolist():
            occurrences = set(row[1:])
            if occurrences != {0}:
                occurrence_df_filtered_list.append(row)

        ## create a dataframe with relative values
        occurrence_df_filtered_relative = pd.DataFrame(occurrence_df_filtered_list)
        occurrence_df_filtered_relative.columns = occurrence_df.columns.tolist()

        ## create a dataframe with absolute values
        occurrence_df_filtered_absolute = pd.DataFrame(occurrence_df_filtered_list)
        occurrence_df_filtered_absolute.columns = occurrence_df.columns.tolist()

        ## convert dataframe to relative occurrence abundance
        for taxon in taxa:
            df = occurrence_df_filtered_relative[["Country", taxon]]
            df_2 = df[taxon]/df[taxon].sum()
            df = df.assign(perc=df_2.values * 100)
            df = df.drop([taxon], axis=1)
            df = df.rename(columns={"perc": taxon})
            occurrence_df_filtered_relative[taxon] = df[taxon]

        fig = go.Figure()
        for row in occurrence_df_filtered_relative.values.tolist():
            occurrences = row[1:]
            country = row[0]
            fig.add_trace(go.Bar(x=list(taxa), y=list(occurrences), text=country, name=country, textposition='auto'))
        fig.update_layout(barmode='stack', width=int(width), height=int(height), template=template)
        fig.update_yaxes(title="GBIF occurrence references (%)")

        n_occurrences = []
        for taxon in taxa:
            n_occurrences.append(len([value for value in occurrence_df_filtered_relative[taxon].values.tolist() if value != 0]))

        fig.add_trace(go.Scatter(x=taxa, y=[105]*len(taxa), text=n_occurrences, name="countries", mode="text"))

        answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
        if answer == "Yes":
            fig.show()
        ## define output files
        output_pdf = Path(str(dirName) + "/" + '_'.join(continents_to_check) + ".pdf")
        output_html = Path(str(dirName) + "/" + '_'.join(continents_to_check) + ".html")
        output_xlsx = Path(str(dirName) + "/" + '_'.join(continents_to_check) + ".xlsx")

        ## write to different sheets, one for absolute data, one for relative
        with pd.ExcelWriter(output_xlsx) as writer:
            occurrence_df_filtered_relative.to_excel(writer, sheet_name='relative', index=False)
            occurrence_df_filtered_absolute.to_excel(writer, sheet_name='absolute', index=False)

        ## write figures
        fig.write_image(str(output_pdf))
        fig.write_html(str(output_html))

        from taxontabletools.create_log import ttt_log
        ttt_log("occurrence analysis", "analysis", TaXon_table_xlsx.name, output_pdf.name, "", path_to_outdirs)

def gbif_occurrence_plot(width, height, template, theme, path_to_outdirs):

    import requests_html, json, glob
    import PySimpleGUI as sg
    import pandas as pd
    from pandas import DataFrame
    import numpy as np
    from pathlib import Path
    import plotly.graph_objects as go
    import os

    # slice function for lists to split up lists
    def slices(list, slice):
        for i in range(0, len(list), slice):
            yield list[i : i + slice]

    ## collect all occurrence tables from the directory
    input_files = sorted(glob.glob(str(path_to_outdirs) + "/Occurrence_analysis/*/*.xlsx"))
    input_files_short = []
    for file in input_files:
        input_files_short.append(str(Path(file).parent).split("/")[-1] + "/" + Path(file).name)

    ##################################################
    # start a second window to ask for the read tables to process

    win2_active = True
    input_files_list = list(slices([sg.Radio(name, "occurrence_tables", default=True) for name in sorted(input_files_short)], 2))
    layout2 = [[sg.Text("Occurrence plots", size=(20,1))],
    [sg.Frame(layout = input_files_list, title = 'Check occurrence table to re-plot')],
    [sg.Button('Plot')],
    [sg.Button('Back')]]

    win2 = sg.Window('Occurrence plots', layout2, keep_on_top=True)
    while True:
        event2, values2 = win2.Read()

        if event2 is None or event2 == 'Back':
            win2.Close()
            win2_active = False
            break

        if event2 == 'Plot':

            if values2 == {}:
                sg.Popup("Please first download occurrence data!")

            else:
                file = Path(str(path_to_outdirs) + "/Occurrence_analysis/" + [input_files_short[value] for key, value in values2.items() if key == True][0])
                occurrence_df_filtered = pd.read_excel(file, sheet_name="relative").fillna(0)
                taxa = occurrence_df_filtered.columns.tolist()[1:]

                fig = go.Figure()
                for row in occurrence_df_filtered.values.tolist():
                    occurrences = row[1:]
                    country = row[0]
                    fig.add_trace(go.Bar(x=list(taxa), y=list(occurrences), text=country, name=country, textposition='auto'))
                fig.update_layout(barmode='stack', width=int(width), height=int(height), template=template)
                fig.update_yaxes(title="GBIF occurrence references (%)")

                n_occurrences = []
                for taxon in taxa:
                    n_occurrences.append(len([value for value in occurrence_df_filtered[taxon].values.tolist() if value != 0]))

                fig.add_trace(go.Scatter(x=taxa, y=[105]*len(taxa), text=n_occurrences, name="countries", mode="text"))

                answer = sg.PopupYesNo('Show plot?', keep_on_top=True)
                if answer == "Yes":
                    fig.show()
                ## define output files
                output_pdf = Path(str(file).replace(".xlsx", ".pdf"))
                output_html = Path(str(file).replace(".xlsx", ".html"))
                fig.write_image(str(output_pdf))
                fig.write_html(str(output_html))

                closing_text = "Occurrence plots are found under:\n" + '/'.join(str(output_pdf).split("/")[-4:])
                sg.Popup(closing_text, title="Finished", keep_on_top=True)

                win2.Close()
                win2_active = False
                break
