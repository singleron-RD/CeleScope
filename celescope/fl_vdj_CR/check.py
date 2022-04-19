"""
from celescope.tools import utils
from celescope.tools.step import Step, s_common
import numpy as np
import pandas as pd
# from bs4 import BeautifulSoup

'''
    # def check(self, sample):
    #     step = 'check'
    #     cmd_line = self.get_cmd_line(step, sample)
    #     rep = f'{self.outdir_dic[sample]["assemble"]}/stat.txt'
    #     rep_bc = f'{self.outdir_dic[sample]["barcode"]}/stat.txt'
    #     html = f'{sample}/{sample}_report.html'
    #     cmd = (
    #         f'{cmd_line} '
    #         f'--rep {rep} '
    #         f'--rep_bc {rep_bc} '
    #         f'--html {html}'
    #     )
    #     self.process_cmd(cmd, step, sample, m=5, x=1)
    
The warning metrics of samples' analysis result are marked in red in the report
'''

warnings = {
    'Q30 of Barcodes':0.90,
    'Q30 of UMIs':0.90,
    'Estimated Number of Cells':10.0,
    'Reads Mapped to Any V(D)J Gene':0.60,
    'Cells With Productive V-J Spanning Pair':0.30,
    'Median used TRA UMIs per Cell':0.1,
    'Median used TRB UMIs per Cell':0.1,
    'Median used IGH UMIs per Cell':0.1,
    'Median used IGL UMIs per Cell':0.1,
    'Median used IGK UMIs per Cell':0.1  
           }

warnings_help = {
    'Q30 of Barcodes': "Ideal > 90%. Application performance may be affected.",
    'Q30 of UMIs': "Ideal > 90%. Application performance may be affected.",
    'Estimated Number of Cells': "Ideal >= 10. This usually indicates poor cell quality, poor library quality, or poor sequencing quality. Application performance is likely to be affected.",
    'Reads Mapped to Any V(D)J Gene': "Ideal > 60%. This can indicate poor specificity of the V(D)J enrichment, use of the wrong germline reference, or the use of an unsupported chemistry type (e.g., using Single Cell 3\' for V(D)J assembly). Application performance may be affected.",
    'Cells With Productive V-J Spanning Pair': "This can indicate poor cell quality, low yield from the RT reaction, poor specificity of the V(D)J enrichment, poor sequencing quality, or the use of an unsupported chemistry type (e.g., using Single Cell 3\' for V(D)J assembly). Application performance may be affected.",
    'Median used TRA UMIs per Cell': "Ideal > 0. This can indicate cells with extremely low TRA expression, poor cell quality, low yield from the RT reaction, or the use of an unsupported chemistry type (e.g., using Single Cell 3\' for V(D)J assembly). Application performance may be affected.",
    'Median used TRB UMIs per Cell': "Ideal > 0. This can indicate cells with extremely low TRB expression, poor cell quality, low yield from the RT reaction, or the use of an unsupported chemistry type (e.g., using Single Cell 3\' for V(D)J assembly). Application performance may be affected.",
    'Median used IGH UMIs per Cell': "Ideal > 0. This can indicate cells with extremely low IGH expression, poor cell quality, low yield from the RT reaction, or the use of an unsupported chemistry type (e.g., using Single Cell 3\' for V(D)J assembly). Application performance may be affected.",
    'Median used IGL UMIs per Cell': "Ideal > 0. This can indicate cells with extremely low IGL expression, poor cell quality, low yield from the RT reaction, or the use of an unsupported chemistry type (e.g., using Single Cell 3\' for V(D)J assembly). Application performance may be affected.",
    'Median used IGK UMIs per Cell': "Ideal > 0. This can indicate cells with extremely low IGK expression, poor cell quality, low yield from the RT reaction, or the use of an unsupported chemistry type (e.g., using Single Cell 3\' for V(D)J assembly). Application performance may be affected."
           }


class Check(Step):
    
    def __init__(self, args, step_name):
        Step.__init__(self, args, step_name)
    
    # input
        self.rep_bc = args.rep_bc
        self.rep = args.rep
        self.seqtype = args.seqtype
        self.html = args.html
        
        self.warning_list=[]


        if self.seqtype =='TCR':
            self.target = [
                'Q30 of Barcodes','Q30 of UMIs','Estimated Number of Cells','Reads Mapped to Any V(D)J Gene','Cells With Productive V-J Spanning Pair',
                'Median used TRA UMIs per Cell','Median used TRB UMIs per Cell']
        elif self.seqtype == 'BCR':
            self.target = [
                'Q30 of Barcodes','Q30 of UMIs','Estimated Number of Cells','Reads Mapped to Any V(D)J Gene','Cells With Productive V-J Spanning Pair',
                'Median used IGH UMIs per Cell','Median used IGL UMIs per Cell','Median used IGK UMIs per Cell']


    @utils.add_log
    def run_check(self):
        rep = dict()
        warning_summary = []

        with open(self.rep_bc,'r') as f:
            lines = (line.strip().split(':') for line in f)
            for line in lines:
                rep[line[0]] = line[1]
        with open(self.rep,'r') as f:
            lines = (line.strip().split(':') for line in f)
            for line in lines:
                rep[line[0]] = line[1]
        
        rep1 = {key:value for key,value in rep.items() if key in self.target}

        for i in rep1.keys():
            rep1[i] = rep1[i].strip()
            if ',' in rep1[i]:
                rep1[i] = rep1[i].replace(',','')
            if '(' in rep1[i]:
                rep1[i] = rep1[i].split('(')[1].strip(')')
            if '%' in rep1[i]:
                rep1[i] = round(float(rep1[i].strip('%'))/100,2)
            else:
                rep1[i] = float(rep1[i])

        for key in rep1:
            if rep1[key] <= warnings[key]:
                warning_summary.append({
                    'item': key,
                    'count': rep1[key],
                    'total_count': np.nan
                })
                
                self.warning_list.append(key)

        if len(self.warning_list) > 0:
            stat_file = self.outdir + '/stat.txt'
            sum_df = pd.DataFrame(warning_summary, columns=['item', 'count', 'total_count'])
            utils.gen_stat(sum_df, stat_file)

    @utils.add_log
    def check_html(self):

        with open(self.html,'r',encoding='utf-8') as file:
            pcontent = file.read()
        sp = BeautifulSoup(pcontent,'html.parser')

        for i in sp.find_all('td'):
            if i.string in self.warning_list:
                new_tag = sp.new_tag("font", color = 'red')
                new_tag.string = i.string
                i.string.replace_with(new_tag)
                new_tag1 = sp.new_tag("font", color = 'red')
                new_tag1.string = i.next_sibling.next_sibling.string
                i.next_sibling.next_sibling.string.replace_with(new_tag1)

        for i in sp.find_all('p'):
            try:
                if i.b.string.strip().strip(':') in self.warning_list:
                    new_tag = sp.new_tag("font", color = 'red')
                    new_tag.string = i.b.string
                    i.b.string.replace_with(new_tag)
                    new_tag1 = sp.new_tag("font", color = 'black')
                    new_tag1.string =  warnings_help[i.b.string.strip().strip(':')]
                    i.b.next_sibling.string.replace_with(new_tag1)
            except StopIteration:
                pass
        
        with open(f'{self.outdir}/{self.sample}_report.html','w',encoding='utf-8') as fp:
            fp.write(sp.prettify())

    def run(self):    
        self.run_check()
        if len(self.warning_list)>0:
            self.check_html()

def check(args):
    step_name = 'check'
    check_obj = Check(args, step_name)
    check_obj.run()

def get_opts_check(parser, sub_program):
    parser.add_argument('--seqtype', help='TCR or BCR', choices=['TCR', 'BCR'], required=True)

    if sub_program:
        parser = s_common(parser)
        parser.add_argument('--rep', help = 'assemble report file', required = True)
        parser.add_argument('--rep_bc', help = 'barcode report file', required = True)
        parser.add_argument('--html', help = 'html report file', required = True)

"""