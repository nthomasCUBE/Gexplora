#
#   Gexplora v0.1d
#
#   2020/03/14: working on the relative genes tagged per total genes in Bin
#   2020/03/15: moving through chromosomes when n_chromosomes>5
#   2020/03/27: integration of stringdb and OMA integration
#   2020/03/28: alternative gene id search for stringdb search
#   2020/03/30: export functions - elements per chromosome
#   2020/04/01: XlsxWriter adaption
#   2020/04/03: adding use case using human annotation file
#   2020/04/06: export functions - elements per bin
#   2020/04/08: fixed ordering of the chromosomes and added title in interfaces
#   2020/04/09: dynamic way to use Bins (e.g. 100, 200, 500, 1000)
#   2020/04/12: adding barley use case
#   2020/04/12: fixing visualization when using different amounts of Bins per chr
#   2020/04/12: adding line charts
#   2020/04/13: if candidate genes are provided, the relative amount per bin is shown [0,100]
#   2020/04/15: adding legend, improved relative frequency visualization
#   2020/04/16: introducing random id when no ID exists in an entry
#   2020/04/26: beginn cis-element analysis panel
#   2020/04/27: begin extracting cis-elements for plants from place
#   2020/10/31: visualization type, reordering menu items
#   2020/10/31: export as image
#   2020/11/14: get sequence and gene trees from gene ids

from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq
import datetime
from datetime import date
from tkinter import *
from tkinter.filedialog import askopenfilename
import tkinter.ttk as ttk
import numpy as np
import requests, sys
import json
import io
import xlsxwriter
from webbrowser import open_new_tab
import urllib
from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont

#
#   Export - Density of elements per Bin
#
def elements_per_bin():
    print("INFO\telements per bin")
    try:
        dens={}
        today = date.today()
        d1 = today.strftime("%Y%m%d")
        fh=open(master.gtf_file, encoding="latin-1")
        workbook = xlsxwriter.Workbook('%s_elements_per_bin.xlsx' % d1)
        worksheet = workbook.add_worksheet()
        worksheet.write("A1","Chromosome")
        worksheet.write("B1","Bin")
        worksheet.write("C1","Amount of elements")
        i1=2
        for dens_ in master.DENS:
            el=master.DENS[dens_]
            c_i=0
            for el_ in el:
                worksheet.write("A"+str(i1),dens_)
                worksheet.write("B"+str(i1),c_i)
                worksheet.write("C"+str(i1),el_)
                c_i=c_i+1
                i1=i1+1
        workbook.close()
    except Exception:
        print("INFO\tgtf file not provided")

#
#   Export - Elements per Chromosome
#
def elements_per_chromosome():
    print("INFO\telements per chromosome")
    try:
        dens={}
        today = date.today()
        d1 = today.strftime("%Y%m%d")
        fh=open(master.gtf_file, encoding="latin-1")
        workbook = xlsxwriter.Workbook('%s_elements_per_chromosome.xlsx' % d1)
        worksheet = workbook.add_worksheet()
        worksheet.write("A1","Chromosome")
        worksheet.write("B1","Amount of elements")
        i1=2
        for line in fh.readlines():
            line=line.strip()
            vals=line.split("\t")
            if(len(vals)>8):
                if(vals[2]==master.tkvar.get()):
                    if(dens.get(vals[0])==None):
                        dens[vals[0]]={}
                    dens[vals[0]][vals[3]]=1
        for dens_ in dens:
            worksheet.write("A"+str(i1),dens_)
            worksheet.write("B"+str(i1),len(dens[dens_].keys()))
            i1=i1+1
        workbook.close()
    except Exception:
        print("INFO\tgtf file not provided")
#
#   StringDB
#
def stringdb_ws():
    print("stringdb_ws")

    my_ids=master.display41.get("1.0","end-1c").split()
    for my_id in my_ids:
        server="http://string-db.org/api/tsv/"
        ext = "abstractsList?identifiers=%s" % my_id
        try:
            r = requests.get(server+ext, headers={ "Content-Type" : "text/tab-separated-values"})
            if not r.ok:
              r.raise_for_status()
              master.display42.insert(END,"SEQUENCE WAS NOT FOUND\n")
              sys.exit()
            else:
                if(len(r.text.split())==1):
                    server="https://rest.ensembl.org/"
                    ext="xrefs/id/%s?content-type=application/json" % my_id
                    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
                    obj=json.loads(r.text)
                    my_id=obj[0]["primary_id"]
                    master.display42.insert(END,"We used following gene name:%s\n" % my_id)
                    server="http://string-db.org/api/tsv/"
                    ext = "abstractsList?identifiers=%s" % my_id
                    r = requests.get(server+ext, headers={ "Content-Type" : "text/tab-separated-values"})
                master.display42.insert(END,r.text)
        except:
            master.display42.insert(END,"UNABLE TO STRINGDB\n")

#
#   OMA
#
def oma_ws():
    from omadb import Client
    c = Client()
    #prot_id = 'P53_RAT'
    my_ids=master.display31.get("1.0","end-1c").split()
    if(len(my_ids)==0):
        master.display32.insert(END,"No orthologous proteins found")
    else:        
        for my_id in my_ids:
            try:
                r = c.proteins[my_id]
                orth=r.orthologs
                master.display32.insert(END,r.orthologs)
            except Exception:
                master.display32.insert(END,"UNABLE TO FIND PROTEIN")
#
#   ENSEMBL
#
def gene_family_ws():
    
    print("INFO\tgene_info_ws")
    server = "https://rest.ensembl.org"
    #my_ids=master.display21.get("1.0","end-1c").split()
    #my_ids_u={}
    #for my_id in my_ids:
    #    my_ids_u[my_id]=1
    my_ids_u=[master.display21.get()]
    print(master.display21.get())
    for my_id in my_ids_u:
        ext = "/family/member/id/%s?content-type=application/json" % my_id
        try:
            print(server+ext)
            r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})
            if not r.ok:
              r.raise_for_status()
              master.display22.insert(END,"SEQUENCE WAS NOT FOUND")
              sys.exit()
            else:
                obj = json.loads(r.text)
                if(len(obj)==0):
                    master.display22.insert(END,"NO GENE FAMILY FOUND")
                else:
                    mem=obj["1"]["members"]
                    for mem_ in mem:
                        psi="-"
                        gen="-"
                        if(mem_.get("protein_stable_id")!=None):
                            psi=mem_["protein_stable_id"]
                        if(mem_.get("genome")!=None):
                            gen=mem_["genome"]
                        if(mem_.get("description")!=None):
                            desc=mem_["description"]
                        master.display22.insert(END,psi+"\t"+gen+"\t"+desc)
                        master.display22.insert(END,"\n")
        except Error as err:
            master.display22.insert(END,err)
    
def gene_info_ws():
    
    print("INFO\tgene_info_ws")
    server = "https://rest.ensembl.org"
    #my_ids=master.display.get()#("1.0","end-1c").split()
    #my_ids_u={}
    #for my_id in my_ids:
    #    my_ids_u[my_id]=1
    my_ids_u=[master.display.get()]
    #print(master.var_1.get())
    for my_id in my_ids_u:
        ext = "/sequence/id/%s?content-type=text/plain" % my_id
        print(ext)
        try: 
            r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})
            if not r.ok:
              r.raise_for_status()
              master.display2.insert(END,"SEQUENCE WAS NOT FOUND")
              sys.exit()
            else:
                master.display2.insert(END,r.text)
        except Error as err:
            master.display2.insert(END,err)

def cis_element_analysis():
    print("cis_element_analysis")

    now = datetime.datetime.today().strftime("%Y%m%d-%H%M%S")
    filename = "cis_element_analysis"+ '.html'
    f = open(filename,'w')

    f2=open('cur_seq.fasta','w')
    my_ids=master.display51.get("1.0","end-1c").split()
    for my_id in my_ids:
        f2.write(my_id+"\n")
    f2.close()

    body="<title>Cis element analysis</title><h1>Results</h1>"

    body=body+"""
    <table border=1>
    <thead>
    <tr>
    <th>Sequence-Id</th>
    <th>Sequence-Match-Position</th>
    <th>PLACE Cis element</th>
    </tr>
    </thead>
    """
    for record in SeqIO.parse('cur_seq.fasta', "fasta"):
        print("%s %i" % (record.id, len(record)))

    wrapper = """<html>
    <head>
    </head>
    <body>%s</body>
    </html>"""

    url="https://www.dna.affrc.go.jp/PLACE/place_dat.shtml";
    file = urllib. request. urlopen(url)
    seqs={}

    for line in file:
        line=str(line).strip()
        if(len(line.split())>2 and line.find("unknown")!=-1):
            vals=line.split()
            if(len(vals)>3):
                seqs[vals[1]]=vals[2]

                for record in SeqIO.parse('cur_seq.fasta', "fasta"):
                    try:
                        instances=[]
                        instances.append(Seq(vals[2],IUPAC.ambiguous_dna))
                        m=motifs.create(instances,alphabet=IUPAC.ambiguous_dna)
                        for pos,seq in m.instances.search(Seq(str(record.seq),alphabet=IUPAC.ambiguous_dna)):
                            body=body+"<tr><td>%s</td><td>%s</td><td>%s</td></tr>" % (record.id,pos,seq);
                    except Exception:
                        pass
    body=body+"<table>"
                    
    url="";
    whole = wrapper % (body)
    f.write(whole)
    f.close()

    #filename = 'file:///Users/username/Desktop/programming-historian/' + filename

    open_new_tab(filename)


    
canvas_width = 600
canvas_height = 600

master = Tk()
master.gois=["---"]
master.gois2=["---"]
master.title("Gexplora - visualizing density of genetic elements in a chromosome");
master.configure(background='lightgreen')

newwin = Toplevel(master, height=10, width=25)
#master.display = Text(newwin, height=5, width=50, bg="lightyellow")
working_list = master.gois #["Option 1", "Option 2", "Option 3", "Option 4"]
master.var_1 = StringVar()
master.display=ttk.Combobox(newwin)#,master.var_1, *working_list)
master.display2 = Text(newwin, height=5, width=50, bg="lightyellow")
button = Button(newwin, text="Get Sequence",command=gene_info_ws)
newwin.destroy()

newwin2 = Toplevel(master, height=10, width=25)
#working_list=master.gois
master.var_1_2=StringVar()
master.display21=ttk.Combobox(newwin2)#,master.var_1_2, *working_list)
#master.display21 = Text(newwin2, height=5, width=50, bg="lightyellow")
master.display22 = Text(newwin2, height=5, width=50, bg="lightyellow")
button2 = Button(newwin2, text="Get Gene Family",command=gene_family_ws)    
newwin2.destroy()

newwin3 = Toplevel(master, height=10, width=25)
master.display31 = Text(newwin3, height=5, width=50, bg="lightyellow")
master.display32 = Text(newwin3, height=10, width=50, bg="lightyellow")
button3 = Button(newwin3, text="Get Orthologous Group in OMA",command=oma_ws)    
newwin3.destroy()

newwin4 = Toplevel(master, height=10, width=25)
master.display41 = Text(newwin4, height=5, width=50, bg="lightyellow")
master.display42 = Text(newwin4, height=10, width=50, bg="lightyellow")
button4 = Button(newwin4, text="StringDB",command=stringdb_ws)    
newwin4.destroy()

newwin5 = Toplevel(master, height=10, width=25)
master.display51 = Text(newwin5, height=5, width=50, bg="lightyellow")
#master.display52 = Text(newwin5, height=10, width=50, bg="lightyellow")
button5 = Button(newwin5, text="Cis element analyis (beta)",command=cis_element_analysis)    
newwin5.destroy()

MAP={}
SYN={}

master.DENS={}

master.CHR_START=0
master.CHR_END=5

def search_genes():
    pass

master.ALL_CHRS={}

#
#   Computes the distribution of gene conservation along the chromosome
#
def calc_distribution(c_type):
    print("INFO\tcalc_distribution\tstart\t%s" % c_type)
    genomes={}
    c_max=0
    gg={}

    if(master.gtf_file==None):
        return;
    
    fh=open(master.gtf_file)
    ri=0
    for line in fh.readlines():
        line=line.strip()
        vals=line.split("\t")
        if(len(vals)==9):
            if(gg.get(vals[0])==None):
                gg[vals[0]]={}
            if(master.ALL_CHRS.get(vals[0])==None):
                master.ALL_CHRS[vals[0]]={}
            master.ALL_CHRS[vals[0]][min(int(vals[3]),int(vals[4]))]=1
            if(len(vals)>0 and vals[2]==c_type):
                genomes[vals[0]]=1
                c_v=min(int(vals[3]),int(vals[4]))
                if(vals[8].find("ID=")!=-1):
                    c_id=vals[8].split("ID=")[1].split(";")[0]
                else:
                    c_id="random_%i" % (ri)
                    ri=ri+1
                gg[vals[0]][c_v]=1
                MAP[c_id]=[vals[0],c_v]
                if(c_v>c_max):
                    c_max=c_v
    print("INFO\tcalc_distribution\tended\t%s" % c_type)

    c_max=-1
    for gg_ in master.ALL_CHRS:
        c_max=max(c_max,max(master.ALL_CHRS[gg_]))
    master.c_step=int(c_max/int(master.tkvar6.get()))

    return(gg)

def OptionMenu_SelectionEvent(event):
    gg=calc_distribution(master.tkvar.get())
    do_calc(gg)

#
#   Genes of interest that should be shown
#
def OptionMenu_SelectionEvent2():    
    filename = askopenfilename()
    fh=open(filename)
    qq={}
    gg=calc_distribution(master.tkvar.get())
    for line in fh.readlines():
        line=line.strip()
        vals=line.split("\t")
        if(len(vals)>0):
            if(MAP.get(vals[0])!=None):
                if(qq.get(MAP[vals[0]][0])==None):
                    qq[MAP[vals[0]][0]]={}
                qq[MAP[vals[0]][0]][MAP[vals[0]][1]]=1
    do_calc(gg, qq)

def cisAnalysis():

    print("INFO\tstringdb")

    newwin5 = Toplevel(master, height=10, width=25)
    newwin5.title("Cis Analysis of Promoters")
    master.display51 = Text(newwin5, height=5, width=50, bg="lightyellow")
    #master.display52 = Text(newwin5, height=10, width=50, bg="lightyellow")
    button = Button(newwin5, text="Analysis",command=cis_element_analysis)

    master.display51.pack() 
    #master.display52.pack() 
    button.pack()

    
def stringdb():

    print("INFO\tstringdb")

    newwin4 = Toplevel(master, height=10, width=25)
    newwin4.title("Obtain Literature referencing the gene")
    master.display41 = Text(newwin4, height=5, width=50, bg="lightyellow")
    master.display42 = Text(newwin4, height=10, width=50, bg="lightyellow")
    button = Button(newwin4, text="Get Sequence",command=stringdb_ws)

    master.display41.pack() 
    master.display42.pack() 
    button.pack()

def gene_info():

    print("INFO\tgene_info")

    newwin = Toplevel(master, height=10, width=25)
    newwin.title("Get Sequence from sequence identifier")
    #master.display = Text(newwin, height=5, width=50, bg="lightyellow")
    working_list = master.gois #["Option 1", "Option 2", "Option 3", "Option 4"]
    master.var_1 = StringVar()
    master.display=ttk.Combobox(newwin)#,master.var_1, *working_list)
    master.display['values']=working_list
    master.display2 = Text(newwin, height=10, width=50, bg="lightyellow")
    button = Button(newwin, text="Get Sequence",command=gene_info_ws)

    master.display.pack() 
    master.display2.pack() 
    button.pack()


def oma():

    print("INFO\tgene_family")

    newwin3 = Toplevel(master, height=10, width=25)
    newwin3.title("Get Orthologous Group")
    master.display31 = Text(newwin3, height=5, width=50, bg="lightyellow")
    master.display32 = Text(newwin3, height=10, width=50, bg="lightyellow")
    button3 = Button(newwin3, text="Orthologous Group",command=oma_ws)    

    master.display31.pack() 
    master.display32.pack() 
    button3.pack()
    
def gene_family():

    print("INFO\tgene_family")

    newwin2 = Toplevel(master, height=10, width=25)
    newwin2.title("Get Gene Family from Gene Identifier")
    working_list=master.gois2
    master.var_1_2 = StringVar()
    master.display21=ttk.Combobox(newwin2)#,master.var_1_2, *working_list)
    #master.display21 = Text(newwin2, height=5, width=50, bg="lightyellow")
    master.display21['values']=working_list
    master.display22 = Text(newwin2, height=10, width=50, bg="lightyellow")
    button2 = Button(newwin2, text="Get Gene Family",command=gene_family_ws)    

    master.display21.pack() 
    master.display22.pack() 
    button2.pack()

def add_chr():
    master.CHR_START=master.CHR_START+5
    master.CHR_END=master.CHR_END+5
    print("INFO\tadd_chr\t%i\t%i" % (master.CHR_START,master.CHR_END))
    gg=calc_distribution(master.tkvar.get())
    do_calc(gg)

def minus_chr():
    master.CHR_START=max(0,master.CHR_START-5)
    master.CHR_END=max(0,master.CHR_END-5)
    if(master.CHR_END==0 and master.CHR_START==0):
        master.CHR_START=0
        master.CHR_END=5
    print("INFO\tminus_chr\t%i\t%i" % (master.CHR_START,master.CHR_END))
    gg=calc_distribution(master.tkvar.get())
    do_calc(gg)

def do_recalc(event):
    print("INFO\tdo_recalc")
    gg=calc_distribution(master.tkvar.get())
    do_calc(gg)

#gg=calc_distribution("one_item")
#c_max=0
#for gg_ in ALL_CHRS:
#    c_max=max(c_max,max(ALL_CHRS[gg_]))
#c_step=int(c_max/BINS)

w = Canvas(master, 
           width=canvas_width, 
           height=canvas_height, bg='lightyellow')

width = 1000
height = 1000
center = height//2
white = (255, 255, 255)
green = (0,128,0)
master.image1 = Image.new("RGB", (width, height), white)
master.w2 = ImageDraw.Draw(master.image1)
w.create_line([0, center, width, center], fill='green')
master.w2.line([0, center, width, center], green)

#master.w2.text((100,100),"medmekl")#,font=font)

def export_canvas():
    print("export_canvas")
    today = date.today()
    d1 = today.strftime("%Y%m%d")
    filename = "%s_canvas.jpg" % d1
    master.image1.save(filename)

export_canvas()


master.tkvar = StringVar(master)
choices = {'---'}
master.tkvar.set('---')
l1=Label(master,text="Element type:",width=30)
master.popupMenu = OptionMenu(master, master.tkvar, *choices, command=OptionMenu_SelectionEvent)

l2=Label(master,text="Candidate genes:",width=30)
popupMenu2=Button(master,text="Find genes (from a text-file)", command=OptionMenu_SelectionEvent2)

master.tkvar3 = StringVar(master)
choices = { 100,90,80,70,60,50,40,30,20,10}
master.tkvar3.set(100)
l4=Label(master,text="Max threshold:",width=30)
popupMenu4 = OptionMenu(master, master.tkvar3, *choices, command=OptionMenu_SelectionEvent)

#tkvar3 = StringVar(master)
#choices = { 100,90,80,70,60,50,40,30,20,10}
#tkvar3.set(100)
l4=Label(master,text="Max threshold:",width=30)
popupMenu51 = Button(master,text="<<", command=minus_chr)
popupMenu52 = Button(master,text=">>", command=add_chr)

master.tkvar6 = StringVar(master)
choices = { 100,150,200,250,300, 500, 1000}
master.tkvar6.set(200)
l6=Label(master,text="Bins total:",width=30)
popupMenu6 = OptionMenu(master, master.tkvar6, *choices, command=do_recalc)

master.var1 = IntVar()
master.var2 = IntVar()
master.var1.set(1)
master.var2.set(1)

def visTypePressed():
    print("visTypePressed")
    master.tkvar71.configure(state=DISABLED)
    master.tkvar72.configure(state=DISABLED)
    gg=calc_distribution(master.tkvar.get())
    do_calc(gg)
    master.tkvar71.configure(state=NORMAL)
    master.tkvar72.configure(state=NORMAL)
    
def get_gtf_file():
    master.gtf_file = askopenfilename()

    master.gois=[]
    master.gois2=[]
    
    choices={}
    fh=open(master.gtf_file, encoding="latin-1")
    for line in fh.readlines():
        line=line.strip()
        vals=line.split("\t")
        if(len(vals)==9):
            choices[vals[2]]=1

            if(vals[2]=="transcript"):
                c_id=vals[8].split("transcript_id ")[1].split()[0]
                if(not(c_id in master.gois)):
                    master.gois.append(c_id.replace("\"","").replace(";",""))
            elif(vals[2]=="mRNA"):
                c_id=vals[8].split("transcript:")[1].split(";")[0]
                if(not(c_id in master.gois)):
                    master.gois.append(c_id.replace("\"","").replace(";",""))
            elif(vals[2]=="gene"):
                try:
                    c_id=vals[8].split("gene_id ")[1].split()[0]
                    if(not(c_id in master.gois2)):
                        master.gois2.append(c_id.replace("\"","").replace(";",""))
                except Exception:
                    try:
                        c_id=vals[8].split("gene:")[1].split(";")[0]
                        if(not(c_id in master.gois2)):
                            master.gois2.append(c_id.replace("\"","").replace(";",""))
                    except:
                        pass
                    
            
    choices=list(choices)
    choices.insert(0,"---")
    master.tkvar = StringVar(master)
    master.tkvar.set('---')
    l1=Label(master,text="Element type:",width=30)
    master.popupMenu = OptionMenu(master, master.tkvar, *choices, command=OptionMenu_SelectionEvent)
    master.popupMenu.grid(row=0,column=1,padx=15)
    
menubar = Menu(master)
filemenu = Menu(master, tearoff=0)
filemenu.add_command(label="Open", command=get_gtf_file)
#filemenu.add_command(label="Save")
filemenu.add_separator()
filemenu.add_command(label="Exit", command=master.quit)
menubar.add_cascade(label="File", menu=filemenu)

filemenu2 = Menu(master, tearoff=0)
filemenu2.add_command(label="Gene Information", command=gene_info)
filemenu2.add_command(label="Gene Family",      command=gene_family)
filemenu2.add_separator()
menubar.add_cascade(label="Ensembl", menu=filemenu2)

filemenu3 = Menu(master, tearoff=0)
filemenu3.add_command(label="Gene Information", command=oma)
filemenu3.add_separator()
menubar.add_cascade(label="OMA", menu=filemenu3)

filemenu4 = Menu(master, tearoff=0)
filemenu4.add_command(label="Obtain PPI", command=stringdb)
filemenu4.add_separator()
menubar.add_cascade(label="StringDB", menu=filemenu4)

filemenu4 = Menu(master, tearoff=0)
filemenu4.add_command(label="Elements per chromosome", command=elements_per_chromosome)
filemenu4.add_command(label="Elements per bin", command=elements_per_bin)
filemenu4.add_command(label="Export to PNG", command=export_canvas)
filemenu4.add_separator()
menubar.add_cascade(label="Export", menu=filemenu4)

filemenu5 = Menu(master, tearoff=0)
filemenu5.add_command(label="Cis-element (input: promoter fasta sequences into textbox)", command=cisAnalysis)
filemenu5.add_separator()
menubar.add_cascade(label="Cis-elements", menu=filemenu5)

master.config(menu=menubar)

l1.grid(row=0,column=0,padx=15)
l2.grid(row=1,column=0,padx=15)
l4.grid(row=3,column=0,padx=15)
l6.grid(row=4,column=0,padx=15)
#master.popupMenu.grid(row=0,column=1,padx=15)
popupMenu2.grid(row=1,column=1,padx=15)
popupMenu4.grid(row=3,column=1,padx=15)
popupMenu6.grid(row=4,column=1,padx=15)
popupMenu51.grid(row=5,column=1,padx=15)
popupMenu52.grid(row=5,column=2,padx=15)
master.tkvar71=Checkbutton(master, text="Linechart", variable=master.var1, command=visTypePressed);
master.tkvar71.grid(row=6, column=1,sticky=W)
master.tkvar72=Checkbutton(master, text="Heatmap", variable=master.var2, command=visTypePressed);
master.tkvar72.grid(row=7, column=1, sticky=W)
w.grid(row=8,column=0,columnspan=2)

def do_calc(gg, qq=None):
    print("INFO\tdo_calc\tstart")
    chrs=[100,80,60,40,30]
    gg_k=gg.keys()

    DENS={}

    if(master.gtf_file==None):
        return;

    # reset of the chromosome view
    y_i=0
    for y in list(master.ALL_CHRS)[0:5]:
        w.create_rectangle(100,100*y_i+10,500,100*y_i+90,fill="white")
        master.w2.line([100,100*y_i+10,500,100*y_i+90],fill="white")
        y_i=y_i+1

    w.delete("all")

    gg_k_sort=sorted(list(gg_k))

    max_nmb_elements=0;

    #
    # getting maximum amount of elements in any bin
    #
    for y in gg_k_sort: #[master.CHR_START:master.CHR_END]:
        for x in range(0,int(master.tkvar6.get())):
            s1=master.c_step*x
            s2=master.c_step*(x+1)
            all_e=gg[y].keys()
            if(len(all_e)>0 and s1<max(all_e)):
                xarr=list(all_e)
                xarr=np.array(xarr)
                xarr=xarr[xarr>=s1]
                xarr=xarr[xarr<=s2]
                cnt=len(xarr)
                if(cnt>max_nmb_elements):
                    max_nmb_elements=cnt
    
    #
    # drawing the headmap on top
    #
    y_i=0
    for y in gg_k_sort[master.CHR_START:master.CHR_END]:
        w.create_rectangle(100,100*y_i+10,300+int(master.tkvar6.get()),100*y_i+90,fill="white")
        master.w2.rectangle([100,100*y_i+10,300+int(master.tkvar6.get()),100*y_i+90],fill="white")
        i=w.create_text(200,100*y_i+85,text="0")
        master.w2.text((200,100*y_i+85),"0",(0,0,0))

        master.DENS[y]=[]

        #print("max_nmb_elements=%s" % max_nmb_elements)

        #
        #   Plotting the bins
        #
        for x in range(0,int(master.tkvar6.get())):
            s1=master.c_step*x
            s2=master.c_step*(x+1)
            all_e=gg[y].keys()
            if(len(all_e)>0 and s1<max(all_e)):
                if(x==0):
                    w.create_text(50,100*y_i+40,text=y)
                    master.w2.text((50,100*y_i+40),str(y),(0,0,0))
                    w.create_text( (300+(int(master.tkvar6.get())-100),100*y_i+85),text=str(round(max(master.ALL_CHRS[y].keys())/1000000,2))+" Mbp")
                    master.w2.text((300+(int(master.tkvar6.get())-100),100*y_i+85),str(round(max(master.ALL_CHRS[y].keys())/1000000,2))+" Mbp",(0,0,0))

                xarr=list(all_e)
                xarr=np.array(xarr)
                xarr=xarr[xarr>=s1]
                xarr=xarr[xarr<=s2]
                cnt=len(xarr)

                cnt2=0
                if(qq!=None):
                    if(qq.get(y)!=None):
                        cnt2=np.array(list(qq[y].keys()))
                        cnt2=len(cnt2[np.isin(cnt2,xarr)])
                        cnt2=(100.0*float(cnt2))/float(cnt)
                else:
                    cnt2=(100.0*float(cnt))/float(max_nmb_elements)
                    
                master.DENS[y].append(cnt2)

                cur_col="white"
                if(cnt2>0.8*float(master.tkvar3.get())):
                    cur_col="red"
                elif(cnt2>0.3*float(master.tkvar3.get())):
                    cur_col="orange"
                elif(cnt2>0.2*float(master.tkvar3.get())):
                    cur_col="green"
                elif(cnt2>0.1*float(master.tkvar3.get())):
                    cur_col="yellow"
                elif(cnt2>0*float(master.tkvar3.get())):
                    cur_col="blue"
                elif(cnt2==0):
                    cur_col="#e0e0d1"
                if(master.var2.get()):
                    w.create_line(200+x,100*y_i+20,200+x,100*y_i+80,fill=cur_col)
                    master.w2.line([200+x,100*y_i+20,200+x,100*y_i+80],fill=cur_col)
                    cur_col="black"

                if(x==0):
                    w.create_rectangle(520,100,540,110,fill="blue")
                    w.create_rectangle(520,90,540,100,fill="yellow")
                    w.create_rectangle(520,80,540,90,fill="green")
                    w.create_rectangle(520,70,540,80,fill="orange")
                    w.create_rectangle(520,60,540,70,fill="red")
                    w.create_text(570,65,text=">80%")
                    w.create_text(570,75,text="[30%,80%]")
                    w.create_text(570,85,text="[20%,30%]")
                    w.create_text(570,95,text="[10%,20%]")
                    w.create_text(570,105,text="]0,10%]")


                    master.w2.rectangle([520,100,540,110],fill="blue")
                    master.w2.rectangle([520,90,540,100],fill="yellow")
                    master.w2.rectangle([520,80,540,90],fill="green")
                    master.w2.rectangle([520,70,540,80],fill="orange")
                    master.w2.rectangle([520,60,540,70],fill="red")
                    master.w2.text((570,65),">80%",(0,0,0))
                    master.w2.text((570,75),"[30%,80%]",(0,0,0))
                    master.w2.text((570,85),"[20%,30%]",(0,0,0))
                    master.w2.text((570,95),"[10%,20%]",(0,0,0))
                    master.w2.text((570,105),"]0,10%]",(0,0,0))

                # not started yet - circular diagrams
                #w.create_arc(100,100,200,200,fill="red")
                #w.create_arc(200,200,100,300,fill="blue")
                #w.create_arc(100,300,0,200,fill="orange")
                #w.create_arc(0,200,100,100,fill="green")

                if(x>0):
                    if(master.var1.get()):
                        w.create_line(200+x-1,100*y_i+100-20-0.6*master.DENS[y][x-1]*0.6,200+x,100*y_i+100-20-0.6*master.DENS[y][x]*0.6,fill=cur_col)                                
                        master.w2.line([200+x-1,100*y_i+100-20-0.6*master.DENS[y][x-1]*0.6,200+x,100*y_i+100-20-0.6*master.DENS[y][x]*0.6],fill=cur_col)                                
            if(x==(int(master.tkvar6.get())-1)):
                if(len(master.DENS[y])>0):
                    w.create_text(350+(int(master.tkvar6.get())-100),100*y_i+25,text=str("max:"+str(round(max(master.DENS[y]),0))))
                    w.create_text(350+(int(master.tkvar6.get())-100),100*y_i+37,text=str("avg:"+str(round(sum(master.DENS[y])/len(master.DENS[y]),0))))
                    w.create_text(350+(int(master.tkvar6.get())-100),100*y_i+49,text=str("min:"+str(round(min(master.DENS[y]),0))))

                    master.w2.text((350+(int(master.tkvar6.get())-100),100*y_i+25),str("max:"+str(round(max(master.DENS[y]),0))),(0,0,0))
                    master.w2.text((350+(int(master.tkvar6.get())-100),100*y_i+37),str("avg:"+str(round(sum(master.DENS[y])/len(master.DENS[y]),0))),(0,0,0))
                    master.w2.text((350+(int(master.tkvar6.get())-100),100*y_i+49),str("min:"+str(round(min(master.DENS[y]),0))),(0,0,0))

        y_i=y_i+1    
    print("INFO\tdo_calc\tended")

master.mainloop()

