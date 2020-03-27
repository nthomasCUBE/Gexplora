#
#   Gexplorr v0.1
#
#   2020/03/14: working on the relative genes tagged per total genes in Bin
#   2020/03/15: moving through chromosomes when n_chromosomes>5
#   2020/03/27: integration of stringdb and OMA integration
#

from tkinter import *
from tkinter.filedialog import askopenfilename
import numpy as np
import requests, sys
import json

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
              master.display42.insert(END,"SEQUENCE WAS NOT FOUND")
              sys.exit()
            else:
                master.display42.insert(END,r.text)
        except:
            master.display42.insert(END,"UNABLE TO STRINGDB")

def oma_ws():
    from omadb import Client
    c = Client()
    #prot_id = 'P53_RAT'
    my_ids=master.display31.get("1.0","end-1c").split()
    for my_id in my_ids:
        r = c.proteins[my_id]
        orth=r.orthologs
        master.display32.insert(END,r.orthologs)

def gene_family_ws():
    
    print("INFO\tgene_info_ws")
    server = "https://rest.ensembl.org"
    my_ids=master.display21.get("1.0","end-1c").split()
    my_ids_u={}
    for my_id in my_ids:
        my_ids_u[my_id]=1
    for my_id in my_ids_u:
        ext = "/family/member/id/%s?content-type=application/json" % my_id
        try: 
            r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})
             
            if not r.ok:
              r.raise_for_status()
              master.display22.insert(END,"SEQUENCE WAS NOT FOUND")
              sys.exit()
            else:
                obj = json.loads(r.text)
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
        except:
            master.display22.insert(END,"UNABLE TO ACCESS ENSEMBL WEBSERVICE\nInternet connection active?")
    
def gene_info_ws():
    
    print("INFO\tgene_info_ws")
    server = "https://rest.ensembl.org"
    my_ids=master.display.get("1.0","end-1c").split()
    my_ids_u={}
    for my_id in my_ids:
        my_ids_u[my_id]=1
    for my_id in my_ids_u:
        ext = "/sequence/id/%s?type=cds" % my_id
        try: 
            r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})
            if not r.ok:
              r.raise_for_status()
              master.display2.insert(END,"SEQUENCE WAS NOT FOUND")
              sys.exit()
            else:
                master.display2.insert(END,r.text)
        except:
            master.display2.insert(END,"UNABLE TO ACCESS ENSEMBL WEBSERVICE\nInternet connection active?")

canvas_width = 600
canvas_height = 600

master = Tk("chromoWIZpy")
master.configure(background='lightgreen')

newwin = Toplevel(master, height=10, width=25)
master.display = Text(newwin, height=5, width=50, bg="lightyellow")
master.display2 = Text(newwin, height=10, width=50, bg="lightyellow")
button = Button(newwin, text="Get Sequence",command=gene_info_ws)
newwin.destroy()

newwin2 = Toplevel(master, height=10, width=25)
master.display21 = Text(newwin2, height=5, width=50, bg="lightyellow")
master.display22 = Text(newwin2, height=10, width=50, bg="lightyellow")
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

BINS=200

MAP={}
SYN={}

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
                c_id=vals[8].split("ID=")[1].split(";")[0]
                gg[vals[0]][c_v]=1
                MAP[c_id]=[vals[0],c_v]
                if(c_v>c_max):
                    c_max=c_v
    print("INFO\tcalc_distribution\tended\t%s" % c_type)

    c_max=-1
    for gg_ in master.ALL_CHRS:
        c_max=max(c_max,max(master.ALL_CHRS[gg_]))
    master.c_step=int(c_max/BINS)

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

def stringdb():

    print("INFO\tstringdb")

    newwin4 = Toplevel(master, height=10, width=25)
    master.display41 = Text(newwin4, height=5, width=50, bg="lightyellow")
    master.display42 = Text(newwin4, height=10, width=50, bg="lightyellow")
    button = Button(newwin4, text="Get Sequence",command=stringdb_ws)

    master.display41.pack() 
    master.display42.pack() 
    button.pack()

def gene_info():

    print("INFO\tgene_info")

    newwin = Toplevel(master, height=10, width=25)
    master.display = Text(newwin, height=5, width=50, bg="lightyellow")
    master.display2 = Text(newwin, height=10, width=50, bg="lightyellow")
    button = Button(newwin, text="Get Sequence",command=gene_info_ws)

    master.display.pack() 
    master.display2.pack() 
    button.pack()


def oma():

    print("INFO\tgene_family")

    newwin3 = Toplevel(master, height=10, width=25)
    master.display31 = Text(newwin3, height=5, width=50, bg="lightyellow")
    master.display32 = Text(newwin3, height=10, width=50, bg="lightyellow")
    button3 = Button(newwin3, text="Orthologous Group",command=oma_ws)    

    master.display31.pack() 
    master.display32.pack() 
    button3.pack()
    
def gene_family():

    print("INFO\tgene_family")

    newwin2 = Toplevel(master, height=10, width=25)
    master.display21 = Text(newwin2, height=5, width=50, bg="lightyellow")
    master.display22 = Text(newwin2, height=10, width=50, bg="lightyellow")
    button2 = Button(newwin2, text="Get Gene Family",command=gene_family_ws)    

    master.display21.pack() 
    master.display22.pack() 
    button2.pack()

def add_chr():
    master.CHR_START=master.CHR_START+5
    master.CHR_END=master.CHR_END+5
    print("INFO\tadd_chr\t%i\t%i" % (master.CHR_START,master.CHR_END))
    gg=calc_distribution(tkvar.get())
    do_calc(gg)

def minus_chr():
    master.CHR_START=max(0,master.CHR_START-5)
    master.CHR_END=max(0,master.CHR_END-5)
    if(master.CHR_END==0 and master.CHR_START==0):
        master.CHR_START=0
        master.CHR_END=5
    print("INFO\tminus_chr\t%i\t%i" % (master.CHR_START,master.CHR_END))
    gg=calc_distribution(tkvar.get())
    do_calc(gg)
    

#gg=calc_distribution("one_item")
#c_max=0
#for gg_ in ALL_CHRS:
#    c_max=max(c_max,max(ALL_CHRS[gg_]))
#c_step=int(c_max/BINS)

w = Canvas(master, 
           width=canvas_width, 
           height=canvas_height, bg='lightyellow')

master.tkvar = StringVar(master)
choices = {'---'}
master.tkvar.set('---')
l1=Label(master,text="Element type:",width=30)
master.popupMenu = OptionMenu(master, master.tkvar, *choices, command=OptionMenu_SelectionEvent)

l2=Label(master,text="Candidate genes:",width=30)
popupMenu2=Button(master,text="Find genes (from a text-file)", command=OptionMenu_SelectionEvent2)

tkvar3 = StringVar(master)
choices = { 100,90,80,70,60,50,40,30,20,10}
tkvar3.set(100)
l4=Label(master,text="Max threshold:",width=30)
popupMenu4 = OptionMenu(master, tkvar3, *choices, command=OptionMenu_SelectionEvent)

tkvar3 = StringVar(master)
choices = { 100,90,80,70,60,50,40,30,20,10}
tkvar3.set(100)
l4=Label(master,text="Max threshold:",width=30)
popupMenu51 = Button(master,text="<<", command=minus_chr)
popupMenu52 = Button(master,text=">>", command=add_chr)

def get_gtf_file():
    master.gtf_file = askopenfilename()

    choices={}
    fh=open(master.gtf_file)
    for line in fh.readlines():
        line=line.strip()
        vals=line.split("\t")
        if(len(vals)==9):
            choices[vals[2]]=1
    master.tkvar = StringVar(master)
    master.tkvar.set('mRNA')
    l1=Label(master,text="Element type:",width=30)
    master.popupMenu = OptionMenu(master, master.tkvar, *choices, command=OptionMenu_SelectionEvent)
    master.popupMenu.grid(row=0,column=1,padx=15)
    
menubar = Menu(master)
filemenu = Menu(master, tearoff=0)
filemenu.add_command(label="Open", command=get_gtf_file)
filemenu.add_command(label="Save")
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

master.config(menu=menubar)

l1.grid(row=0,column=0,padx=15)
l2.grid(row=1,column=0,padx=15)
l4.grid(row=3,column=0,padx=15)
#master.popupMenu.grid(row=0,column=1,padx=15)
popupMenu2.grid(row=1,column=1,padx=15)
popupMenu4.grid(row=3,column=1,padx=15)
popupMenu51.grid(row=4,column=1,padx=15)
popupMenu52.grid(row=4,column=2,padx=15)
w.grid(row=5,column=0,columnspan=2)

def do_calc(gg, qq=None):
    print("INFO\tdo_calc\tstart")
    chrs=[100,80,60,40,30]
    gg_k=gg.keys()

    if(master.gtf_file==None):
        return;

    # reset of the chromosome view
    y_i=0
    for y in list(master.ALL_CHRS)[0:5]:
        w.create_rectangle(100,100*y_i+10,500,100*y_i+90,fill="white")
        y_i=y_i+1

    w.delete("all")
    
    # drawing the headmap on top
    y_i=0
    for y in list(gg_k)[master.CHR_START:master.CHR_END]:
        w.create_rectangle(100,100*y_i+10,500,100*y_i+90,fill="white")
        i=w.create_text(200,100*y_i+85,text="0")
        
        for x in range(0,BINS):
            s1=master.c_step*x
            s2=master.c_step*(x+1)
            all_e=gg[y].keys()
            if(len(all_e)>0 and s1<max(all_e)): # and max(ALL_CHRS[y])>15000000):
                if(x==0):
                    w.create_text(50,100*y_i+40,text=y)
                    w.create_text(400,100*y_i+85,text=str(round(max(master.ALL_CHRS[y].keys())/1000000,2))+" Mbp")

                xarr=list(all_e)
                xarr=np.array(xarr)
                xarr=xarr[xarr>=s1]
                xarr=xarr[xarr<=s2]
                cnt=len(xarr)

                cnt2=0
                if(qq!=None):
                    if(qq.get(y)!=None):
                        cnt2=np.array(list(qq[y].keys()))
                        if(cnt>0):
                            cnt2_N=len(cnt2[np.isin(cnt2,xarr)])
                            cnt2=round(100.0*cnt2_N/len(xarr))
                    cnt=cnt2_N
                else:
                    cnt2=0

                cur_col="white"
                if(cnt>0.8*float(tkvar3.get())):
                    cur_col="red"
                elif(cnt>0.3*float(tkvar3.get())):
                    cur_col="orange"
                elif(cnt>0.2*float(tkvar3.get())):
                    cur_col="green"
                elif(cnt>0.1*float(tkvar3.get())):
                    cur_col="yellow"
                elif(cnt>0*float(tkvar3.get())):
                    cur_col="blue"
                elif(cnt==0):
                    cur_col="#e0e0d1"
                w.create_line(200+x,100*y_i+20,200+x,100*y_i+80,fill=cur_col)
        y_i=y_i+1    
    print("INFO\tdo_calc\tended")

#do_calc(gg)

master.mainloop()

