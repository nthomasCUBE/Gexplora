#
#   2002/03/10
#

from tkinter import *
from tkinter.filedialog import askopenfilename

import numpy as np

import requests, sys

canvas_width = 600
canvas_height = 600

master = Tk()

BINS=200

MAP={}
SYN={}

def search_genes():
    pass

ALL_CHRS={}

#
#   Computes the distribution of gene conservation along the chromosome
#
def calc_distribution(c_type):
    print("INFO\tcalc_distribution\tstart\t%s" % c_type)
    genomes={}
    c_max=0
    gg={}
    fh=open("brachypodium1.2_20100223_MIPSGFF.gff")
    for line in fh.readlines():
        line=line.strip()
        vals=line.split("\t")
        if(gg.get(vals[0])==None):
            gg[vals[0]]={}
        if(ALL_CHRS.get(vals[0])==None):
            ALL_CHRS[vals[0]]={}
        ALL_CHRS[vals[0]][int(vals[3])]=1
        if(len(vals)>0 and vals[2]==c_type):
            genomes[vals[0]]=1
            c_v=max(int(vals[3]),int(vals[4]))
            c_id=vals[8].split("ID=")[1].split(";")[0]
            gg[vals[0]][c_v]=1
            MAP[c_id]=[vals[0],c_v]
            if(c_v>c_max):
                c_max=c_v
    print("INFO\tcalc_distribution\tended\t%s" % c_type)
    return(gg)

def OptionMenu_SelectionEvent(event):
    gg=calc_distribution(tkvar.get())
    do_calc(gg)

#
#   Genes of interest that should be shown
#
def OptionMenu_SelectionEvent2():
    filename = askopenfilename()
    fh=open(filename)
    gg={}
    for line in fh.readlines():
        line=line.strip()
        vals=line.split("\t")
        if(len(vals)>0):
            if(MAP.get(vals[0])!=None):
                if(gg.get(MAP[vals[0]][0])==None):
                    gg[MAP[vals[0]][0]]={}
                gg[MAP[vals[0]][0]][MAP[vals[0]][1]]=1
    do_calc(gg)


#
#   Determines whether a gene has an ortholog or not
#
def OptionMenu_SelectionEvent3():
    filename = askopenfilename()
    fh=open(filename)
    SYN={}
    for line in fh.readlines():
        line=line.strip()
        vals=line.split("\t")
        if(len(vals)>0):
            SYN[vals[0]]=vals[1]


def gene_info():
    print("INFO\tgene_info")
    server = "https://rest.ensembl.org"
    ext = "/sequence/id/BRADI_5g16797v3?type=cds"
     
    r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})
     
    if not r.ok:
      r.raise_for_status()
      sys.exit()
     
    print(r.text)     

gg=calc_distribution("one_item")
c_max=0
for gg_ in ALL_CHRS:
    c_max=max(c_max,max(ALL_CHRS[gg_]))
c_step=int(c_max/BINS)

w = Canvas(master, 
           width=canvas_width, 
           height=canvas_height)

tkvar = StringVar(master)
choices = { 'mRNA','exon','one_item'}
tkvar.set('one_item')
popupMenu = OptionMenu(master, tkvar, *choices, command=OptionMenu_SelectionEvent)
popupMenu.pack()

popupMenu2=Button(master,text="Find genes (from a text-file)", command=OptionMenu_SelectionEvent2);
popupMenu2.pack()

popupMenu3=Button(master,text="orthology (from a text-file)", command=OptionMenu_SelectionEvent3);
popupMenu3.pack()

tkvar3 = StringVar(master)
choices = { 100,90,80,70,60,50,40,30,20,10}
tkvar3.set(100)
popupMenu3 = OptionMenu(master, tkvar3, *choices, command=OptionMenu_SelectionEvent)
popupMenu3.pack()

menubar = Menu(master)
filemenu = Menu(master, tearoff=0)
filemenu.add_command(label="Open")
filemenu.add_command(label="Save")
filemenu.add_separator()
filemenu.add_command(label="Exit", command=master.quit)
menubar.add_cascade(label="File", menu=filemenu)

filemenu2 = Menu(master, tearoff=0)
filemenu2.add_command(label="Gene Information", command=gene_info)
filemenu2.add_separator()
menubar.add_cascade(label="Ensembl", menu=filemenu2)

master.config(menu=menubar)

w.pack()

def do_calc(gg):
    print("INFO\tdo_calc\tstart")
    chrs=[100,80,60,40,30]
    gg_k=gg.keys()

    # reset of the chromosome view
    y_i=0
    for y in ALL_CHRS:
        w.create_rectangle(100,100*y_i+10,500,100*y_i+90,fill="white")
        y_i=y_i+1
        
    # drawing the headmap on top
    y_i=0
    for y in gg_k:
        w.create_rectangle(100,100*y_i+10,500,100*y_i+90,fill="white")
        w.create_text(200,100*y_i+85,text="0")
        for x in range(0,BINS):
            s1=c_step*x
            s2=c_step*(x+1)
            all_e=gg[y].keys()
            if(len(all_e)>0 and s1<max(all_e) and max(ALL_CHRS[y])>15000000):
                if(x==0):
                    w.create_text(50,100*y_i+40,text=y)
                    w.create_text(400,100*y_i+85,text=str(round(max(ALL_CHRS[y].keys())/1000000,2))+" Mbp")
                xarr=list(all_e)
                xarr=np.array(xarr)
                xarr=xarr[xarr>=s1]
                xarr=xarr[xarr<=s2]
                cnt=len(xarr)

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
                    cur_col="black"
                w.create_line(200+x,100*y_i+20,200+x,100*y_i+80,fill=cur_col)
        y_i=y_i+1    
    print("INFO\tdo_calc\tended")

do_calc(gg)

mainloop()

