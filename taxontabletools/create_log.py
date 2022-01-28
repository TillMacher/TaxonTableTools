import os
import pandas as pd
import numpy as np
from datetime import datetime
from pathlib import Path
import networkx as nx
import matplotlib.pyplot as plt
import PySimpleGUI as sg

def ttt_log(task, type, input_file, output_file, parameters, path_to_outdirs):

    dateTimeObj = datetime.now()
    time = str(dateTimeObj)

    log_xlsx = Path(str(path_to_outdirs) + "/log.xlsx")

    if os.path.exists(Path(log_xlsx)):
        log_df = pd.read_excel(log_xlsx)
        log_df = pd.DataFrame(log_df.values.tolist() + [[task, input_file, output_file, parameters, time, type]],columns=log_df.columns.tolist())
    else:
        log_df = pd.DataFrame([[task, input_file, output_file, parameters, time, type]],columns=["task", "input", "output", "parameters", "time", "type"])

    log_df.to_excel(log_xlsx, index=False)

def ttt_log_network(path_to_outdirs):

    log_xlsx = Path(str(path_to_outdirs) + "/log.xlsx")

    if not os.path.exists(Path(log_xlsx)):
        sg.PopupError("No log file exists!", title="Error", keep_on_top=True)
        raise RuntimeError

    win2_layout = [
                [sg.Text('',size=(1,1))],
                [sg.Text('Please adjust the plot size to the size of the network.')],
                [sg.Text("Plot size (w,h):"), sg.Input("12", size=(3,1), key="x_log"), sg.Input("12", size=(3,1), key="y_log")],
                [sg.Text('',size=(1,1))],
                [sg.Button('Continue', key='Continue')],
                ]

    win2 = sg.Window('Projects', win2_layout)
    event, values = win2.read()
    x_log = int(values["x_log"])
    y_log = int(values["y_log"])
    win2.close()

    df = pd.read_excel(log_xlsx)

    edge_label_dict = {}
    node_label_dict = {}
    for task in df.values.tolist():
        # collect data from the log table
        input = task[1]
        output = task[2]
        edge_label = task[0]
        type = task[-1]
        # create a dict entry for the edge name
        edge_label_dict[(input, output)] = edge_label
        # createa a dict entry for the node name
        if type == "processing":
            node_label_dict[input] = input
            node_label_dict[output] = output
        else:
            node_label_dict[input] = input
            node_label_dict[output] = ""

    plt.figure(figsize=(x_log, y_log))
    G=nx.from_pandas_edgelist(df, 'input', 'output', create_using=nx.Graph())
    nx.draw(G, with_labels=False, node_color="skyblue", node_size=50, width=3, edge_color = "skyblue", edge_cmap=plt.cm.Blues, pos=nx.kamada_kawai_layout(G), font_size=8)
    nx.draw_networkx_edge_labels(G, pos=nx.kamada_kawai_layout(G), edge_labels=edge_label_dict, font_color='blue', font_size=6)
    nx.draw_networkx_labels(G, pos=nx.kamada_kawai_layout(G), labels=node_label_dict, font_color='black', font_size=8)
    plt.axis('off')
    plt.show(block=False)

    answer = sg.PopupYesNo('Save figure?', keep_on_top=True)
    if answer == "Yes":
        log_pdf = Path(str(path_to_outdirs) + "/log.pdf")
        plt.savefig(log_pdf)
        plt.close()
        closing_text = "Log network is found under:\n" + log_pdf.name
        sg.Popup(closing_text, title="Finished", keep_on_top=True)
    else:
        plt.close()
