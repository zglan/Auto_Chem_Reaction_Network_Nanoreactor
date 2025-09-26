import argparse

class ArgParse:
    
    @staticmethod
    def get_args():
        parser = argparse.ArgumentParser()
        
        ### xyz format trajectory file ###
        parser.add_argument(
            '-i', '--inputfile', 
            help='Input trajectory file')

        ### reaction event dectection paraments ###
        parser.add_argument(
            '--interval', 
            help='Set interval', type=float, default=1)
        parser.add_argument(
            '--scale',
            help='Set frame in detect reaction, default 30 frames', 
            type=int, default=30)
        parser.add_argument(
            '--mode',
            help='Select mode: reac or prod, default extract from product', 
            default='prod')
        
        ### HMM paraments ###
        parser.add_argument(
            '--hidmatr', 
            help='Matrix hidden of HMM parameters', 
            type=float, nargs=4, default=[0.999, 0.001, 0.001, 0.999])
        parser.add_argument(
            '--obvmatr', 
            help='Matrix obvserve of HMM parameters', 
            type=float, nargs=4, default=[0.6, 0.4, 0.4, 0.6])

        ### Spin and Charge set ###
        parser.add_argument(
            '--tspin', 
            help='Set spin for entire system', 
            type=int, default=1)
        parser.add_argument(
            '--tchrg', 
            help='Set charge for entire system', 
            type=int, default=0)

        ### G16 job set and paraments ###
        parser.add_argument(
            '--jobset',
            help='Set job set, default scale is -30~30, interval is 10', 
            type=int, nargs=2, default=[30, 10])
        parser.add_argument(
            '--ts',
            help='Proceed gaussion ts job', action='store_true')
        parser.add_argument(
            '--opt', 
            help='Proceed gaussion opt job', action='store_true')
        
        ### Orca(NEB) XTB(rsmd-pp, opt) job set ###
        parser.add_argument(
            '--smooth',
            help='Proceed smooth trajectory job to get double end points. \
                !! Need nebterpolator !!', 
            action='store_true')
        parser.add_argument(
            '--neb', 
            help='Proceed orca neb-ts job', action='store_true')
        parser.add_argument(
            '--xtbpath', 
            help='Proceed xtb path job', action='store_true')
        parser.add_argument(
            '--xtbopt', 
            help='Proceed xtb opt job', action='store_true')
        
        ### network analysis ###
        parser.add_argument(
            '--network', 
            help='Proceed network analysis job', action='store_true')    

        ### draw reactions ###
        parser.add_argument(
            '--draw', 
            help='draw reactions', action='store_true')

        ### refine reactions ###
        parser.add_argument(
            '--refine', 
            help='refine reactions', action='store_true') 
        
        ### load reaction pre-analysis ###
        parser.add_argument(
            '--load', 
            help='load reaction pre-analysis', action='store_true')  

        return parser