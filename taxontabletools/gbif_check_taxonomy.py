def gbif_parent_check(phylum_name, taxon_name, taxonomy_check):

    import requests_html, json
    import PySimpleGUI as sg
    import pandas as pd
    from pandas import DataFrame
    import numpy as np
    from pathlib import Path

    #taxon_name = "Radix ampla" taxon_level = "species"
    ## create an html session
    with requests_html.HTMLSession() as session:
        ## generate html request name
        request_name = '%20'.join(taxon_name.split(' '))
        ## request that name
        r = session.get("https://api.gbif.org/v1/species/match?verbose=true&name=" + request_name + "&limit=1")
        ## parse json
        res = json.loads(r.text)

        gbif_result = []

        if 'note' in res.keys():
            if "Multiple equal matches" in res['note']:
                for match in res["alternatives"]:
                    if phylum_name == match["phylum"]:
                        try:
                            for taxon_level in taxonomy_check:
                                    gbif_result.append(match[taxon_level.lower()])
                            return gbif_result
                            break
                        except:
                            gbif_result.append("")
                            return gbif_result
                            break
            else:
                for taxon_level in taxonomy_check:
                    try:
                        gbif_result.append(res[taxon_level.lower()])
                    except:
                        try:
                            gbif_result.append(res["alternatives"][0][taxon_level.lower()])
                        except:
                            gbif_result.append("")
                return gbif_result
        else:
            gbif_result.append("")
            return gbif_result

def gbif_check_taxonomy(TaXon_table_xlsx, path_to_outdirs):

    import requests_html, json
    import PySimpleGUI as sg
    import pandas as pd
    from pandas import DataFrame
    import numpy as np
    from pathlib import Path

    TaXon_table_xlsx = Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx)
    taxon_levels = ["Phylum","Class","Order","Family","Genus","Species"]
    OTUs_list = TaXon_table_df["ID"].values.tolist()

    taxonomy_check_dict = {}

    ############################################################################
    ## create the progress bar window
    layout = [[sg.Text('Progress bar')],
              [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
              [sg.Cancel()]]
    window_progress_bar = sg.Window('Progress bar', layout, keep_on_top=True)
    progress_bar = window_progress_bar['progressbar']
    progress_update = 0
    progress_increase = 1000 / len(OTUs_list) + 1
    ############################################################################

    for OTU in TaXon_table_df[["Phylum","Class","Order","Family","Genus","Species"]].fillna("").values.tolist():

        for i, taxonomy in enumerate(OTU):
            if taxonomy == "":
                phylum_name = OTU[0]
                taxon_name = OTU[i-1]
                taxonomy_check = taxon_levels[0: i]
                result = gbif_parent_check(phylum_name, taxon_name, taxonomy_check)
                query = OTU[0: i]
                if query != result:
                    if len(query) != 6:
                        add = 6 - len(query)
                        query = query + [''] * add
                    if len(result) != 6:
                        add = 6 - len(result)
                        result = result + [''] * add
                    query = ",".join(query)
                    taxonomy_check_dict[query] = result
                break

            elif i == 5:
                phylum_name = OTU[0]
                taxon_name = OTU[5]
                taxonomy_check = taxon_levels
                result = gbif_parent_check(phylum_name, taxon_name, taxonomy_check)
                if OTU != result:
                    query = ",".join(OTU)
                    taxonomy_check_dict[query] = result

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

    TaXon_table_list = []

    for OTU in TaXon_table_df.fillna("").values.tolist():
        taxonomy = OTU[1:7]
        search_key = ','.join(taxonomy)
        if (search_key in taxonomy_check_dict.keys() and taxonomy_check_dict[search_key] != ['']*6):
            replacement_taxonomy = taxonomy_check_dict[search_key]
            replacement_OTU = [OTU[0]] + replacement_taxonomy + OTU[7:]
            TaXon_table_list.append(replacement_OTU)
        else:
            TaXon_table_list.append(OTU)

    file_name = TaXon_table_xlsx.stem
    output_name = Path(str(path_to_outdirs) + "/TaXon_tables/" + file_name + "_gbif" + ".xlsx")
    df_new = pd.DataFrame(TaXon_table_list, columns=(TaXon_table_df.columns.values.tolist()))
    df_new.to_excel(output_name, sheet_name = 'TaXon table', index=False)

    change_log_list = []
    for key, value in taxonomy_check_dict.items():
        change_log_list.append(["Input:"] + key.split(","))
        change_log_list.append(["Gbif:"] + value)

    change_log_df = pd.DataFrame(change_log_list, columns=(["Change"] + taxon_levels))
    change_log_name = Path(str(path_to_outdirs) + "/GBIF/" + file_name + "_gbif_log" + ".xlsx")
    change_log_df = pd.DataFrame(change_log_list, columns=(["Change"] + taxon_levels))
    change_log_df.to_excel(change_log_name, sheet_name = 'TaXon table', index=False)

    closing_text = "Taxon table is found under:\n" + '/'.join(str(output_name).split("/")[-4:]) + "\n\n" + "Log file is found under:\n" + '/'.join(str(change_log_name).split("/")[-4:])
    sg.Popup(closing_text, title="Finished", keep_on_top=True)

    from taxontabletools.create_log import ttt_log
    ttt_log("gbif check", "processing", TaXon_table_xlsx.name, output_name.name, "nan", path_to_outdirs)
