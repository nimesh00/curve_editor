'''
                        MIT License

    Copyright (c) 2020 Nimesh Khandelwal

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.


                    ME751A FINAL PROJECT
                    Submitted by: Nimesh Khandelwal
                    Team Size: 1
                    Instructor: Prof. Anupam Saxena

*   Software Requirements: Python3.5+(Developed and tested on 3.8.6), numpy, tkinter(python3-tk: Developed and tested on version (3.8.6-1))

*   To install the required softwares, run the following comands:
    sudo apt install python3
    sudo apt install python3-numpy
    sudo apt install python3-tk
            OR
    # you can use pip3 as well to install numpy and tk packages. A simple google search will yeild exact commands for the same.

*   To run the program, run: python3 curveEditor_2D.py # make sure you have installed the required libraries.

*   This is a basic curve editor program that is capable of calculating and
    dynamically manipulating Beizer and B-spline segments for given control
    points.

*   Avoid Using higher degree as the program tends to slow down considerably for degree greater than 5.

*   Red Curve is the Spline Curve and Blue curve is the Bezeir Curve

*   Functions of the button(s) present in the program:
    ADD:    Add the control point at the coordinates specified by the X and Y values given in the entry fields above this button.
            You can only add integers values in the given Range: X = [0, 900] Y = [0, 400]
    UPDATE: Update the coordinates of the control points selected from the list below this button. To select a control point just click on that.
            The control point cooridnates get updated to values in the X and Y entry field above.
    DELETE: Delete the Selected control point. To select the control point you need to click on it in the list shown.

    (You can scroll in the list to view all control points. If no control point is selected, using UPDATE and DELETE act on the top-most control point in the list.)

    SHOW/HIDE BEZEIR CURVE: Show/Hide Bezeir Curve calculated via de Casteljau scheme.
                            This currently has a bug that to hide it the first time you have to click twice. Works correctly afterwards.
    SHOW/HIDE SPLINE CURVE: Show/Hide B-Spline Curve calculated via DeBoor scheme.
                            This currently has a bug that to hide it the first time you have to click twice. Works correctly afterwards.
    RESET: Deletes all the control points and clears canvas.

* References:
    1) Computer Aided Engineering Design - Anupam Saxena, Birendra Sahay; Springer 2005
    2) Curves and Surfaces In Geometric Modeling: Theory And Algorithms - Gallier, Jean, 2018
    3) https://docs.python.org/3/library/tkinter.html
    4) https://github.com/giuluck/CurvesEditor
    5) https://nbviewer.jupyter.org/github/empet/geom_modeling/blob/master/FP-Bezier-Bspline.ipynb (For FP implementation of deBoor algorithm)
    6) https://stackoverflow.com/
    7) https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node18.html


Author(s): Nimesh Khandelwal
e-mail: nimesh6798@gmail.com
GitHub: https://github.com/nimesh00

'''

from tkinter import *
from tkinter import ttk
import numpy as np
from functools import partial

# FP based implementation of DeBoor algorithm
cvx=lambda a, b, t: (1-t)*a+t*b

def cvxP(P_Q_t):
    P, Q, t=P_Q_t
    return (cvx(P[0], Q[0], t), cvx(P[1], Q[1], t))


def cvxCtrlP(b, t):
    return list(map(cvxP, zip(b[:-1], b[1:], [t]*(len(b)-1))))


def omega(u, k, t):
    return map(lambda j: (t-u[j]) / (u[j+k]-u[j]),  range(1,k+1))

def cvxList(d, alpha):
    return  list(map(cvxP, zip(d[:-1], d[1:], alpha)))

def DeBoor(d, u, k, t):
    if(len(d)==1):
        return d[0]
    else:
        return  DeBoor(cvxList(d, omega(u, k,t)), u[1:-1], k-1,t)

def Bspline(d, k=3, N=100, type=2):
    L=len(d)-1
    n=L+k+1
    #extend the control polygon
    if type == 0:
        d = [d[0] for j in range (k)] + d
    elif type == 1:
        d += [d[-1] for j in range (k)]
    elif type == 2:
        dupIndices = int((k + 1) / 2)
        d = [d[0] for j in range (dupIndices)] + d + [d[-1] for j in range (dupIndices)]
    else:
        d += [d[j] for j in range (k+1)]

    #make uniform knots
    u=np.arange(n+k+1)
    # define the period T
    T=u[n]-u[k]
    #extend the sequence of knots
    u[:k]=u[n-k:n]-T
    u[n+1:n+k+1]=u[k+1:2*k+1]+T
    u=list(u)

    curve=[]
    for J in range(k, n+1):

        if u[J]<u[J+1]:
           t=np.linspace(u[J], u[J+1], N)
        else: continue

        curve+=[DeBoor(d[J-k:J+1],  u[J-k:J+k+1], k, t[j])  for j in range (N) ]
    return curve


# Imperetive implementation of de Casteljau scheme
def deCasteljauImp(b,t):
    N=len(b)
    if(N<2):
        raise InvalidInputError("The  control polygon must have at least two points")
    a=np.copy(b)
    for r in range(1,N):
        a[:N-r,:]=(1-t)*a[:N-r,:]+t*a[1:N-r+1,:]
    return a[0,:]

def BezierCv(b, nr=200):# compute nr points on the Bezier curve of control points in list bom
    t=np.linspace(0, 1, nr)
    return [deCasteljauImp(b, t[k]) for k in range(nr)]



# FP based Implementation of de Casteljau scheme
def deCasteljauF(b, t):
    # de Casteljau scheme - computes the point p(t) on the Bezier curve of control polygon b
    if len(b)>1:
        return deCasteljauF(cvxCtrlP( b, t), t)
    else:
        return b[0]

def BezCurve(b, nr=200):
    #computes nr points on the Bezier curve of control points b

    return map(lambda s: deCasteljauF(b,s),  map(lambda j: j*1.0/(nr-1), xrange(nr)))

class Curve2D:
    def __init__(self, master):
        self.controlPointLabels = []
        self.controlPoints = []
        self.controlPointInt = []
        self.controlPointID = []
        self.controlLines = []
        self.showPointXLabel = []
        self.showPointYLabel = []
        self.controlPointDeleteButton = []
        self.bezeirCurve = []
        self.splineCurve = []
        self.degree = 3
        # self.bezeirCurvePoints = []
        self.showControlPointRowStart = 2
        self.master = master
        master.title("2D Curve Manipulator")
        self.width = 900
        self.height = 400
        self.canvas = Canvas(master, width=self.width, height=self.height, background="white")
        # self.create_canvas_grid()

        self.addPointLabel = Label(master, text="Enter X and Y coordinates of the point")

        self.addPointEntryXLabel = Label(master, text="X: ")
        self.addPointEntryYLabel = Label(master, text="Y: ")


        vcmdx = master.register(self.validatex)
        vcmdy = master.register(self.validatey)
        self.addPointEntryX = Entry(master, validate="key", validatecommand=(vcmdx, '%P'))
        self.addPointEntryY = Entry(master, validate="key", validatecommand=(vcmdy, '%P'))

        Button(master,text=" Add  ",command=self.addEntry).grid(row = 2, column = 3)
        Button(master,text="Update",command=self.updateEntry).grid(row = 2, column = 7)
        Button(master,text="Delete",command=self.deleteEntry).grid(row = 2, column = 8)


        self._drag_data = {"x": 0, "y": 0, "item": None}

        self.canvas.tag_bind("controlPoint", "<ButtonPress-1>", self.drag_start)
        self.canvas.tag_bind("controlPoint", "<ButtonRelease-1>", self.drag_stop)
        self.canvas.tag_bind("controlPoint", "<B1-Motion>", self.drag)

        self.scroll = Scrollbar(master, orient=VERTICAL)
        self.select = Listbox(master, yscrollcommand=self.scroll.set, height=6)
        self.scroll.config (command=self.select.yview)
        self.scroll.grid(row=3, column=3)
        self.select.grid(row=3, column=3, columnspan = 4)

        vcmdd = master.register(self.validated)
        Button(self.master, text="Show/Hide Bezeir Curve", command=self.toggleBezeir).grid(row = 4, column = 0)
        Button(self.master, text="Show/Hide Spline Curve", command=self.togglespline).grid(row = 5, column = 0)
        Label(self.master, text="Spline Type: ").grid(row = 5, column = 1)
        self.splineType = StringVar()
        self.splineTypeList = ttk.Combobox(master, width = 30, textvariable = self.splineType)
        self.splineTypeList['values'] = ('Open, Left End Clamped', 'Open, Right End Clamped', 'Open, Both Ends Clamped', 'Closed Curve')
        self.splineTypeList.grid(row=5, column = 2)
        self.splineTypeList.current(2)
        self.splineTypeList.bind("<<ComboboxSelected>>", self.splineTypeSelectUpdate)
        Label(self.master, text="Degree: ").grid(row = 5, column = 3)
        self.degreeEntry = Entry(self.master, validate="key", validatecommand=(vcmdd, '%P'))
        self.degreeEntry.grid(row = 5, column = 4)

        Button(master, text="RESET", command = self.reset_everything).grid(row = 5, column = 5)


        # LAYOUT

        self.canvas.grid(row=0, column = 0, rowspan = 4, columnspan = 3)
        self.addPointLabel.grid(row = 0, column = 3, columnspan = 4)
        self.addPointEntryXLabel.grid(row = 1, column = 3, sticky=W)
        self.addPointEntryYLabel.grid(row = 1, column = 5, sticky=W)
        self.addPointEntryX.grid(row = 1, column = 4, sticky=W)
        self.addPointEntryY.grid(row = 1, column = 6, sticky=W)


    def splineTypeSelectUpdate(self, event):
        self.setSelect()

    def calculateSpline(self):
        d = [[float(self.controlPoints[i][0].get()), float(self.controlPoints[i][1].get())] for i in range(len(self.controlPoints))]
        return Bspline(d, k = self.degree, N = 100, type = self.splineTypeList.current())
    def calculateBezeir(self):
        b = [[float(self.controlPoints[i][0].get()), float(self.controlPoints[i][1].get())] for i in range(len(self.controlPoints))]
        return BezierCv(b)

    def togglespline(self):
        if len(self.controlPoints) < 2:
            return

        if len(self.splineCurve) == 0:
            # print("Calculating first time!!")
            self.splineCurve = self.plot_curve(self.calculateSpline(), width = 2, color = 'red')
        else:
            # print("Toggling state!!")
            if self.canvas.itemcget(self.splineCurve[0],'state') == 'normal':
                itemstate = 'hidden'
                self.degreeEntry.grid_remove()
                self.splineTypeList.grid_remove()
            else:
                itemstate = 'normal'
                self.degreeEntry.grid()
                self.splineTypeList.grid()
            for curvePeice in self.splineCurve:
                self.canvas.itemconfigure(curvePeice, state=itemstate)
        self.setSelect()

    def toggleBezeir(self):
        if len(self.controlPoints) < 2:
            return
        if len(self.bezeirCurve) == 0:
            # print("Calculating first time!!")
            self.bezeirCurve = self.plot_curve(self.calculateBezeir(), width = 2, color = 'blue')
        else:
            # print("Toggling state!!")
            if self.canvas.itemcget(self.bezeirCurve[0],'state') == 'normal':
                itemstate = 'hidden'
            else:
                itemstate = 'normal'
            for curvePeice in self.bezeirCurve:
                self.canvas.itemconfigure(curvePeice, state=itemstate)
        self.setSelect()

    def plot_curve(self, p,width = 1, color = 'black'):
        curve = []
        for i in range(len(p)):
            if i < len(p) - 1:
                curve.append(self.canvas.create_line(p[i][0], self.height - p[i][1], p[i + 1][0], self.height - p[i + 1][1], width = width, fill=color))
        return curve

    def reset_everything(self):
        while len(self.controlPoints) > 0:
            self.deleteEntry()

    def create_canvas_grid(self):
        for line in range(0, self.width, 10): # range(start, stop, step)
            self.canvas.create_line([(line, 0), (line, self.height)], fill="#DCDCDC", tags='grid_line_w')
        for line in range(0, self.height, 10):
            self.canvas.create_line([(0, line), (self.width, line)], fill="#DCDCDC", tags='grid_line_h')

    def setSelect (self) :
        self.select.delete(0,END)
        for controlPoint in self.controlPoints:
            self.select.insert (END, str(controlPoint[0].get() + ', ' + controlPoint[1].get()))
        self.update()

    def whichSelected (self) :
        try:
            return int(self.select.curselection()[0])
        except:
            if len(self.controlPoints) != 0:
                return 0;
            else:
                return

    def addEntry(self) :
        control_point_x = StringVar()
        control_point_y = StringVar()
        control_point_x.set(str(self.entered_x))
        control_point_y.set(str(self.entered_y))
        controlPoint = [control_point_x, control_point_y]
        self.controlPoints.append(controlPoint)
        self.plot_points([[self.entered_x, self.entered_y]], radius = 5)
        if len(self.controlPoints) > 1:
            self.controlLines.append(self.canvas.create_line(float(self.controlPoints[-2][0].get()), self.height - float(self.controlPoints[-2][1].get()), \
                                                             float(self.controlPoints[-1][0].get()), self.height - float(self.controlPoints[-1][1].get())))
        self.setSelect ()

    def updateEntry(self) :
        control_point_x = StringVar()
        control_point_y = StringVar()
        control_point_x.set(str(self.entered_x))
        control_point_y.set(str(self.entered_y))
        self.controlPoints[self.whichSelected()] = [control_point_x, control_point_y]
        self.canvas.coords(self.controlPointID[self.whichSelected()], self.entered_x, self.height - self.entered_y, self.entered_x + 5, self.height - self.entered_y + 5)
        self.setSelect ()

    def deleteEntry(self):
        # print(self.whichSelected())
        if len(self.controlPoints) == 0:
            return
        self.canvas.delete(self.controlPointID[self.whichSelected()])
        del self.controlPointID[self.whichSelected()]
        if self.whichSelected() == 0:
            if len(self.controlLines) != 0:
                self.canvas.delete(self.controlLines[self.whichSelected()])
                del self.controlLines[self.whichSelected()]
        elif self.whichSelected() == len(self.controlPoints) - 1:
            self.canvas.delete(self.controlLines[self.whichSelected() - 1])
            del self.controlLines[self.whichSelected() - 1]
        else:
            self.canvas.delete(self.controlLines[self.whichSelected()])
            self.canvas.delete(self.controlLines[self.whichSelected() - 1])
            del self.controlLines[self.whichSelected()]
            del self.controlLines[self.whichSelected() - 1]
            self.controlLines.insert(self.whichSelected(), self.canvas.create_line(float(self.controlPoints[self.whichSelected() - 1][0].get()), self.height - float(self.controlPoints[self.whichSelected() - 1][1].get()), \
                                                                                   float(self.controlPoints[self.whichSelected() + 1][0].get()), self.height - float(self.controlPoints[self.whichSelected() + 1][1].get())))
        del self.controlPoints[self.whichSelected()]
        self.setSelect ()


    def validatex(self, new_text):
        if not new_text: # the field is being cleared
            self.entered_x = 0
            return True

        try:
            self.entered_x = int(new_text)
            return True
        except ValueError:
            return False

    def validatey(self, new_text):
        if not new_text: # the field is being cleared
            self.entered_y = 0
            return True

        try:
            self.entered_y = int(new_text)
            return True
        except ValueError:
            return False
    def validated(self, new_text):
        if not new_text: # the field is being cleared
            self.degree = 0
            return True

        try:
            self.degree = int(new_text)
            return True
        except ValueError:
            return False


    def plot_points(self, points, radius = 1):
        for point in points:
            x1, y1 = point[0], self.height - point[1]
            x2, y2 = x1 + radius, y1 + radius
            self.controlPointID.append(self.canvas.create_oval(x1, y1, x2, y2, fill="#000000", outline="#000000", tags=("controlPoint",)))


    def drag_start(self, event):
        """Begining drag of an object"""
        # record the item and its location
        self._drag_data["item"] = self.canvas.find_closest(event.x, event.y)[0]
        self._drag_data["x"] = event.x
        self._drag_data["y"] = event.y

    def drag_stop(self, event):
        """End drag of an object"""
        # reset the drag information
        self._drag_data["item"] = None
        self._drag_data["x"] = 0
        self._drag_data["y"] = 0

    def drag(self, event):
        """Handle dragging of an object"""
        # compute how much the mouse has moved
        delta_x = event.x - self._drag_data["x"]
        delta_y = event.y - self._drag_data["y"]
        # move the object the appropriate amount
        self.canvas.move(self._drag_data["item"], delta_x, delta_y)
        # record the new position
        self._drag_data["x"] = event.x
        self._drag_data["y"] = event.y
        self.setSelect()

    def update(self):
        # print("control Points, IDs, Lines: ", len(self.controlPoints), len(self.controlPointID), len(self.controlLines))
        # print(self.controlPointID)
        for i in range(len(self.controlPointID)):
            self.controlPoints[i][0].set(str(self.canvas.coords(self.controlPointID[i])[0]))
            self.controlPoints[i][1].set(str(self.height - self.canvas.coords(self.controlPointID[i])[1]))
            if i < len(self.controlPointID) - 1:
                self.canvas.coords(self.controlLines[i], float(self.controlPoints[i][0].get()), self.height - float(self.controlPoints[i][1].get()), \
                                                         float(self.controlPoints[i + 1][0].get()), self.height - float(self.controlPoints[i + 1][1].get()))
        if len(self.controlPoints) < 2:
            for i in range(len(self.bezeirCurve)):
                self.canvas.delete(self.bezeirCurve[i])
            del self.bezeirCurve[:]
            for i in range(len(self.splineCurve)):
                self.canvas.delete(self.splineCurve[i])
            del self.splineCurve[:]
            return

        # Visible curve updation
        if len(self.bezeirCurve) != 0:
            self.bezeirCurvePoints = self.calculateBezeir()
            if len(self.bezeirCurve) != len(self.bezeirCurvePoints) - 1:
                for i in range(len(self.bezeirCurve)):
                    self.canvas.delete(self.bezeirCurve[i])
                del self.bezeirCurve[:]
                self.bezeirCurve = self.plot_curve(self.bezeirCurvePoints, width = 2, color = 'blue')
            else:
                if len(self.bezeirCurvePoints) > 0:
                    for i in range(len(self.bezeirCurve)):
                        self.canvas.coords(self.bezeirCurve[i], self.bezeirCurvePoints[i][0], self.height - self.bezeirCurvePoints[i][1], \
                                                                self.bezeirCurvePoints[i + 1][0], self.height - self.bezeirCurvePoints[i + 1][1])
        if len(self.splineCurve) != 0:
            self.splineCurvePoints = self.calculateSpline()
            if len(self.splineCurve) != len(self.splineCurvePoints) - 1:
                for i in range(len(self.splineCurve)):
                    self.canvas.delete(self.splineCurve[i])
                del self.splineCurve[:]
                self.splineCurve = self.plot_curve(self.splineCurvePoints, width = 2, color = 'red')
            else:
                if len(self.splineCurvePoints) > 0:
                    for i in range(len(self.splineCurve)):
                        self.canvas.coords(self.splineCurve[i], self.splineCurvePoints[i][0], self.height - self.splineCurvePoints[i][1], \
                                                                self.splineCurvePoints[i + 1][0], self.height - self.splineCurvePoints[i + 1][1])



root = Tk()
my_gui = Curve2D(root)
root.mainloop()
