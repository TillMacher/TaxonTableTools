def convert_to_presence_absence(TaXon_table_xlsx, path_to_outdirs):

    from pathlib import Path
    import PySimpleGUI as sg
    import  pandas as pd

    TaXon_table_xlsx =  Path(TaXon_table_xlsx)
    TaXon_table_df = pd.read_excel(TaXon_table_xlsx, header=0)
    # create presence/absence table
    presence_absence_list = []
    for col in TaXon_table_df.values.tolist():
        presence_absence_list.append(col[0:10] + [int(1) if reads != 0 else int(0) for reads in col[10:]])
    df_pa = pd.DataFrame(presence_absence_list)
    df_pa.columns = TaXon_table_df.columns.tolist()
    output_file = Path(str(path_to_outdirs) + "/" + "TaXon_tables" + "/" + TaXon_table_xlsx.stem + "_pa.xlsx")
    df_pa.to_excel(output_file, index=False, sheet_name = 'TaXon table')
    closing_text = "Presence absence tables is found in: " + str(path_to_outdirs) + "/TaXon_tables/"
    sg.Popup(closing_text, title="Finished", keep_on_top=True)

    from taxontabletools.create_log import ttt_log
    ttt_log("presence absence conversion", "processing", TaXon_table_xlsx.name, output_file.name, "nan", path_to_outdirs)
