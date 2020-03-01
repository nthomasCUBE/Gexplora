from tkinter import *

canvas_width = 600
canvas_height = 600

master = Tk()

def calc_distribution():
    genomes={}
    c_max=0
    gg={}
    fh=open("brachypodium1.2_20100223_MIPSGFF.gff")
    for line in fh.readlines():
        line=line.strip()
        vals=line.split("\t")
        if(len(vals)>0):
            genomes[vals[0]]=1
            c_v=max(int(vals[3]),int(vals[4]))
            if(gg.get(vals[0])==None):
                gg[vals[0]]={}
            gg[vals[0]][c_v]=1
            if(c_v>c_max):
                c_max=c_v
    return(gg)

gg=calc_distribution()
c_max=0
for gg_ in gg:
    c_max=max(c_max,max(gg[gg_]))
c_step=int(c_max/200)

w = Canvas(master, 
           width=canvas_width, 
           height=canvas_height)

tkvar = StringVar(master)
choices = { 'All Chromosomes','Chromosome-1','Chromosome-2','Chromosome-3','Chromosome-4','Chromosome-5'}
tkvar.set('All Chromosomes')
popupMenu = OptionMenu(master, tkvar, *choices)
popupMenu.pack()

w.pack()

chrs=[100,80,60,40,30]

gg_k=gg.keys()
y_i=0
for y in gg_k:
    w.create_rectangle(100,100*y_i+10,500,100*y_i+90,fill="white")
    for x in range(0,200):
        print("\t%s\t%i" % (y,x))
        s1=c_step*x
        s2=c_step*(x+1)

        if(s2>max(all_e) and max(all_e)>1000000):
            cnt=0
            all_e=gg[y].keys()
            all_E={}
            for z in range(s1,s2):
                all_E[z]=1
            all_E_arr=all_E.keys()
            all_E_arr=list(set(all_E_arr) & set(all_e))
            cnt=len(all_E_arr)

            cur_col="white"
            if(cnt>100):
                cur_col="red"
            elif(cnt>50):
                cur_col="orange"
            elif(cnt>20):
                cur_col="green"
            elif(cnt>10):
                cur_col="yellow"
            elif(cnt>0):
                cur_col="blue"
            elif(cnt==0):
                cur_col="black"
            w.create_line(200+x,100*y_i+20,200+x,100*y_i+80,fill=cur_col)
    y_i=y_i+1    
mainloop()

