import numpy as np
import toml, impedance, logging, argparse

logging.basicConfig(level=logging.DEBUG)

parser = argparse.ArgumentParser()
parser.add_argument('input')
parser.add_argument('-c','--code')
parser.add_argument('-n','--name')
parser.add_argument('-p','--pretty-print', action='store_true')
parser.add_argument('-C','--print-configuration', action='store_true')
parser.add_argument('-s','--single-line', action='store_true')
parser.add_argument('-l','--list-lines', action='store_true')
parser.add_argument('-P','--python', action='store_true')
parser.add_argument('-j','--julia', action='store_true')


args = parser.parse_args()
configs = impedance.load_configs(args.input)

def print_mat(Z):
    print('[')
    for zi in Z:
        print('  [ ', end='')

        for zj in zi:
            print(f'{zj:9.6f}, ', end='')

        print('],')
    
    print(']')

for key, cfg  in configs['lines'].items():
    calc = {key: cfg}
    cfg = configs['lines'][key]

    if args.list_lines:

        if 'name' in cfg:
            print(f"{key} name = {cfg['name']}")
        else:
            print(f"{key}")
        continue

    if args.code and key != args.code:
        continue

    if 'cable' in cfg:
        cfg['construction'] = 'underground'

        if args.print_configuration:
            print(toml.dumps(calc))

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
        if args.print_configuration:
            print(toml.dumps(cfg))

        Zij, Zin, znn = impedance.overhead_line(cfg)
        Z = impedance.kron(Zij, Zin, znn)
        Zs = np.array2string(Z, max_line_width=80)
        
        if args.pretty_print:
            print(f'Impedance:\n{Zs}')
        elif args.single_line:
            print('[impedance]')
            print(f'R = {np.round(Z.real, decimals=6).tolist()}')
            print(f'X = {np.round(Z.imag, decimals=6).tolist()}')
            print()            
        else:
            print('[impedance]')
            print('R = ', end='')
            print_mat(Z.real)

            print('\nX = ', end='')
            print_mat(Z.imag)
            print()


        try:
            z0, z1 = impedance.sequences(Z)
            if args.pretty_print:
                print(f'\n\nSequence impedances:\nz1 = {z1}\nz0 = {z0}')
            elif args.python:
                print(f'\n\nSequence impedances:')
                print(f'z1 = {z1.real:0.6f} + {z1.imag:0.6f}j')
                print(f'z0 = {z0.real:0.6f} + {z0.imag:0.6f}j')
            elif args.julia:
                print(f'\n\nSequence impedances:')
                print(f'z1 = {z1.real:0.6f} + {z1.imag:0.6f}im')
                print(f'z0 = {z0.real:0.6f} + {z0.imag:0.6f}im')                
            else:
                print('[sequence_impedances]')
                print(f'r1 = {z1.real:0.6f}')
                print(f'x1 = {z1.imag:0.6f}')      
                print(f'r0 = {z0.real:0.6f}')
                print(f'x0 = {z0.imag:0.6f}')                              
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





