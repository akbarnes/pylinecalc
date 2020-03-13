import PySimpleGUI as sg  
import numpy as np
import toml, impedance, logging

logging.basicConfig(level=logging.DEBUG)

configs = impedance.load_configs("line_configs.toml", long_line_keys=True, duplicate_line_keys=False)

line_keys = list(configs['lines'].keys())

sg.ChangeLookAndFeel('SystemDefaultForReal')      


line_list = sg.Listbox(values=line_keys, size=(40, 25), key='lines')
line_col = [[sg.Text('Lines', size=(10, 1))], [line_list], [sg.Button('Calculate'), sg.Button('Reload')]]
cfg_col = [[sg.Text('Config', size=(10,1))], [sg.Multiline(size=(40, 30), key='config')], [sg.Button('Calculate Config'), sg.Button('Save Config')]]
results_col = [[sg.Text('Results', size=(10,1))], [sg.Multiline(size=(60, 30), key='results')], [sg.Button('Save Results')]]

menu_def = [['&File', ['&Open', '&Save', '---', 'Properties', 'E&xit'  ]],
            ['&Edit', ['Paste', ['Special', 'Normal',], 'Undo'],],
            ['&Help', '&About...'],]

layout = [[sg.Menu(menu_def)], [sg.Column(line_col), sg.Column(cfg_col), sg.Column(results_col)]]

window = sg.Window('Imp', layout)      

while True:                             # The Event Loop
    event, values = window.read() 
    print(event, values)       

    if event in (None, 'Exit'):      
        break      

    if event != 'Calculate':
        continue

    if len(values['lines']) == 0:
        continue

    key = values['lines'][0]

    if key not in configs['lines']:
        continue 

    cfg = configs['lines'][key]


    if 'cable' in cfg:
        cfg['construction'] = 'underground'
        cfg_str = toml.dumps(cfg)
        window['config'].update(cfg_str)

        Zpp, Zpn, Znn, Dpp, Dpn, Dnn = impedance.underground_line(cfg)
        Z = impedance.kron(Zpp, Zpn, Znn)
        Zs = np.array2string(Z, max_line_width=80)
        result_str = f'Impedance:\n{Zs}'

        try:
            z0, z1 = impedance.sequences(Z)
            result_str += f'\n\nSequence impedances:\nz1 = {z1}\nz0 = {z0}'
        except:
            logging.error('Error calculating sequence impedances')


        try:
            Y = impedance.underground_line_shunt(cfg)
            Ys = np.array2string(Y, max_line_width=80)
            result_str += f'\n\nShunt Admittance:\n{Y}'
        except:
            logging.erro('Error calculating shunt admittance')

    else:
        cfg['construction'] = 'overhead'
        cfg_str = toml.dumps(cfg)
        window['config'].update(cfg_str)

        Zij, Zin, znn = impedance.overhead_line(cfg)
        Z = impedance.kron(Zij, Zin, znn)
        Zs = np.array2string(Z, max_line_width=80)
        result_str = f'Impedance:\n{Zs}'

        try:
            z0, z1 = impedance.sequences(Z)
            result_str += f'\n\nSequence impedances:\nz1 = {z1}\nz0 = {z0}'
        except:
            logging.error('Error calculating sequence impedances')

        try:
            Pij, Pin, pnn = impedance.overhead_line_shunt(cfg)
            P = impedance.kron(Pij, Pin, pnn)
            C = impedance.capacitance(P)
            Cs = np.array2string(C, max_line_width=80)
            result_str += f'\n\nCapacitance:\n{Cs}'
        except:
            logging.error('Error calculating shunt capacitances')

    window['results'].update(result_str)




window.close()
