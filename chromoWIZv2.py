from tkinter import *

canvas_width = 600
canvas_height = 600

master = Tk()

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

for y in range(0,5):
    w.create_rectangle(100,100*y+10,500,100*y+90,fill="white")
    for x in range(200,200+chrs[y]*2):
        w.create_line(x,100*y+20,x,100*y+80,fill="blue")
    
mainloop()

