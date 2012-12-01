#!/usr/bin/env python

try:
    import SVGdraw
except ImportError:
    pass

"""
Working solely in the 'diagonal' system from here on!
"""
class Point:
    def __init__(self,x,y,knot):
        if not isinstance(x,int) and isinstance(y,int) and not (x+y)%2:
            raise Exception("Point must be ints and add to an even number")
        self.x=x
        self.y=y
        self.knot=knot

    def __eq__(self,other):
        if isinstance(other,Point):
            return self.x==other.x and self.y==other.y and self.knot==other.knot
        # should be instance or 2-element list; exceptions otherwise.
        return self.x==other[0] and self.y==other[1]

    def __ne__(self,other):
        return not self.__eq__(other)

    def __repr__(self):
        # Don't show the knot.  Maybe the modulus.
        return "Point(%d,%d)(%d)"%(self.x, self.y, self.knot.xmodulus)

    def __str__(self):
        return self.__repr__()

    def __hash__(self):
        # I need this for set manipulation.  Is this safe to mess with?
        return hash(self.knot)+1000*self.x+10*self.y+10

    @classmethod
    def Pointlist(cls, knot, lst):
        rv=[]
        for p in lst:
            if isinstance(p,Point):
                rv.append(Point(p.x,p.y,knot))
            else:
                # Let exceptions happen here, if they will.
                rv.append(Point(p[0],p[1],knot))
        return rv

    def lines(self):
        # y might exceed xmodulus.  It does make sense in this case
        # to take a y-coordinate modulo the xmodulus.
        return [(self.x+self.y%self.knot.xmodulus)%self.knot.xmodulus,
                (self.x-self.y%self.knot.xmodulus)%self.knot.xmodulus]

class Knot:
    def __init__(self, ptlist):
        if not ptlist:
            return
        ptlist=Point.Pointlist(self,ptlist)
        self.pivots=ptlist[:]   # do I need to make a copy, you think?
        self.valid=True
        self.reconfigure()

    def reconfigure(self):
        minx=min([p.x for p in self.pivots])
        miny=min([p.y for p in self.pivots])
        if not (minx,miny) in self.pivots:
            raise Exception("Knot must have a lower left corner")
        # Normalize lower left corner to (0,0)
        for i in range(0,len(self.pivots)):
            self.pivots[i].x -= minx
            self.pivots[i].y -= miny

        # sort, row-first I think.
        self.pivots.sort(cmp=(lambda a,b: cmp(a.y,b.y) or cmp(a.x,b.x)))
        # Looks confusing.  modulus needs to be +1 for odd leads, +2 for even.
        self.xmodulus=max([p.x for p in self.pivots])
        self.xmodulus += 2 - self.xmodulus%2
        self.ymax=max([p.y for p in self.pivots]) # useful to know.
        self.validate()

    def validate(self):
        # Warn about these?
        # (a) odd length of self.pivots
        # (b) xmodulus should == len(self.pivots)-1
        # (c) points not on even-sum lattices? (checked in Point)
        # (d) lines that don't connect?
        if len(self.pivots)%2 or self.xmodulus != len(self.pivots):
            self.valid=False
        # Other validations?

    def __repr__(self):
        return "Knot(%s)"%self.pivots.__repr__()
    def __str__(self):
        return self.__repr__()

    def line(self,p,slope):
        # Return all points that share a line with p, *not* including
        # p itself (p should be a pivot--error if not?  Or allow).
        # There really should be only one.
        # slope is the direction, the slope of the line in question: +1 or -1
        
        # Lines are returned as intercepts, in (-,+) order.
        slope=0 if slope <= 0 else 1
        intercept=p.lines()[slope]
        rv=filter(lambda pt: intercept==pt.lines()[slope], self.pivots)
        rv.remove(p)            # error if not there... catch it?
        return rv

    def strands(self):
        "How many strands will it take to tie this knot?"
        visited=set(self.pivots)
        howmany=0
        while visited:
            howmany+=1
            start=None
            p=visited.pop()
            # Put it back, so the test works out okay.
            visited.add(p)
            direction=1
            while p!=start:
                if start is None:
                    start=p
                nxt=self.line(p,direction)
                if len(nxt) != 1:
                    # This knot is inconsistent with tying at all!
                    return 0
                nxt=nxt[0]
                visited.remove(p)
                p=nxt
                direction*=-1
        return howmany

    def circuit(self,p):
        "Return list of pivots making circuit starting at point p"
        if not p in self.pivots:
            raise Exception("Start point must be pivot")
        start=None
        direction=1
        rv=[]
        while p!=start:
            print "Circuit at %s"%p
            if start is None:
                start=p
            rv.append(p)
            nxt=self.line(p,direction)
            if len(nxt) != 1:
                raise Exception("Knot is not tyable.")
            p=nxt[0]
            direction*=-1
        return rv

    def pointsbetween(self,start,end):
        "Return list of lattice-points on line between start and end, *not* including the endpoints."
        slope=self.slopebetween(start,end)
        if not slope:
            # Points do not share a (diagonal) line
            return []
        rv=[]
        # How do I tell which direction to go in?  Is it a crime to
        # just trial-and-error?
        direction=1
        while True:
            (x,y)=((start.x+direction)%self.xmodulus,
                   (start.y+slope*direction))
            # non-end-point can't be at x=0, y=0, or y=ymax, so
            # safe to do modulus there.
            while end != (x,y) and x>0 and y>0 and y<self.ymax:
                rv.append(Point(x,y,self))
                x+=direction    # Slope only affects one coordinate
                x%=self.xmodulus
                y+=slope*direction
            if end != (x,y):
                if direction==1: # It better be
                    direction = -1
                    continue
                else:           # Can't really happen
                    raise Exception("Could not complete line %s--%s"%
                                    (str(start),str(end)))
            else:
                break
        return rv

    def slopebetween(self,p1,p2):
        "Return slope of line between two points; zero if there is no line."
        l1=p1.lines()
        l2=p2.lines()
        if l1[0]==l2[0]:
            return -1
        if l1[1]==l2[1]:
            return 1
        return 0

    def oncircuit(self,path):
        "Returns the lattice points, _in order_, along the given path, _not_ including pivots.  If the path can be closed, we assume it is to be closed."
        rv=[]
        for i in range(0,len(path)-1): # stop after getting TO the end.
            pts=self.pointsbetween(path[i],path[i+1])
            if not pts:
                return []       # nothing for bad paths!
            rv.extend(pts)
        # At the end, though, close if possible.  Otherwise
        # pointsbetween will just return [], which disappears.
        rv.extend(self.pointsbetween(path[-1],path[0]))
        return rv

    def svgout(self,stroke_width=0.3,scale=100,circle_radius=0.3,
               startat=None,segpercolor=None,crossings=True):
        try:
            if type(SVGdraw)!=type(__builtins__):
                return None
        except NameError:
            return None

        cols=['#000000', 
              '#800000', '#808000', '#008080', '#000080',
              '#ff2000', '#ffff20', '#20ffff', '#0020ff',
              '#ff0080', '#ff8000', '#8000ff', '#80ff00']
        colordiv=len(self.pivots/6)
        if segpercolor is not None and segpercolor > 0:
            colordiv=segpercolor
        def coloriter():
            colcounter=0
            while True:
                yield cols[int(colcounter/colordiv)%len(cols)]
                colcounter+=1

        color=coloriter()

        svg=SVGdraw.svg(width=(max_x-min_x+2)*scale,
                        height=(max_y-min_y+2)*scale,
                        x=(min_x-1)*scale,
                        y=(min_y-1)*scale,
                        transform="scale(%d)"%scale)
        defs=SVGdraw.defs(id="defs")
        svg.addElement(defs)
        maingroup=SVGdraw.group(id="main")
        svg.addElement(maingroup)
        # Positive slopes and negative slopes.
        plus=SVGdraw.group(id="plus")
        minus=SVGdraw.group(id="minus")
        circuit=self.circuit(self.pivots[0])
        for i in range(0,len(circuit)-1):
            here=circuit[i]
            nxt=circuit[i+1]
            col=color()
            # Check for wraparound!
            # ???
            # XXXXXXX
            line=SVGdraw.line(here.x,here.y,there.x,there.y,
                                   stroke_width=stroke_width,
                                   stroke=color)
