#!/usr/bin/env python

try:
    import SVGdraw
except ImportError:
    pass

import sys

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

    def lines(self,nowrap=False):
        "Return (y(-1), y(1)), the y-intercepts of the lines passing through this point of slope -1 and 1 respectively."
        # (I guess these are also the negatives of the x-intercepts)
        # It does make sense in this case to take a y-coordinate modulo the
        # xmodulus.
        if (nowrap):
            return [(self.x+self.y),
                    (self.x-self.y)]
        mod=self.knot.xmodulus
        return [(self.y+self.x)%mod,
                (self.y-self.x)%mod]

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

    @classmethod
    def TH(cls,leads,bights):
        "Return a simple LxB Turks-Head"
        return cls(zip(range(0,2*bights,2),[0]*bights) +
                   zip(range((leads%2),2*bights,2),[leads]*bights))

    @classmethod
    def Layers(cls,layers):
        """Try to build a flat-bottomed multi-tier TH.
Knot.Layers(layers):

layers is a list of the form [[n1,h1],[n2,h2],[n3,h3]...]
where n1 is the number of pivots to be placed at hight h1 and so on.
Total length of the knot, number of bights or pivots along the (flat)
bottom is the sum of all the n's.  Each n must divide the total 
(have gcd(n,total)>1).  Trial-and-error attempt to build single-strand
by varying phase?
"""
        from fractions import gcd
        total=sum([e[0] for e in layers])
        if any([gcd(total,e[0])<2 for e in layers]):
            raise Exception("Layers must not be prime to total.")
        pivots=[]
        layertable=[]
        for layer in layers:
            [number,height]=layer
            # if you're filling in 6 out of 18, the row looks like 6
            # sections of length 3, each with one place filled.  If
            # it's 12 out of 18, still 6 sections of length 3, with
            # two places filled.  etc.
            sections=gcd(total,number)
            size=total/sections
            howmany=number/sections
            thisrow=[]
            # Odd rows need to be on odd lattice points.
            shift=height%2
            for i in range(0,total,size):
                for j in range(0,howmany):
                    thisrow.append((2*(i+j)+shift,height))
            layertable.append(thisrow)
        # At this point, should have layers but all with phase=0.
        # Need to make trials now?
        def shiftover(row, shift):
            "Return a row that looks like this one only shifted over, modulo the total."
            return zip([(p[0]+shift)%(2*total) for p in row],
                       [p[1] for p in row])
        import itertools
        def flatten(lol):
            # flatten list of lists (only one level deep)
            return list(itertools.chain(*lol))
        def shifted(ltbl,shifts):
            """Take a layer table ltbl.  Assemble a list by concatenating 
            ltbl[0] shifted over by shifts[0] with
            ltbl[1] shifted over by shifts[1] and
            ltbl[2] shifted over by shifts[2] and so forth.
            Return that.
            """
            ss=[]
            for i in range(0,len(ltbl)):
                ss+=shiftover(ltbl[i],shifts[i])
            return ss
        # Try all the possibilities?  Ow.
        ranges=[xrange(0,i) for i in [len(l) for l in layertable]]
        # LEN(L) FOR L IN LAYERTABLE IS WRONG!!! ^^^^^^^^^^^
        itr=itertools.product(*ranges)
        bottom=zip(range(0,2*total,2),[0]*total)
        results=set()
        for sft in itr:
            k=Knot(bottom+shifted(layertable,sft))
            strnd=0
            try:
                strnd=len(k.strands())
            except Exception:
                pass
            if strnd:
                results.add(k)
        return results
                    

    def __repr__(self):
        return "Knot(%s)"%str([(p.x,p.y) for p in self.pivots])
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

    def strands(self, start=None):
        "Return circuits that cover ALL pivots.  i.e. a three-stranded knot will have a list of three circuits"
        visited=set(self.pivots)
        rv=[]
        while visited:
            if start is None:
                start=visited.pop()
                visited.add(start)
            current=self.circuit(start)
            visited.difference_update(set(current))
            rv.append(current)
            start=None
        return rv

    def circuit(self,p=None):
        "Return list of pivots making circuit starting at point p"
        if p is None:
            p=self.pivots[0]    # you deserve an exception if pivots is []
        if not p in self.pivots:
            raise Exception("Start point must be pivot")
        start=p
        direction=1
        rv=[]
        while p!=start or not rv:
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
        # In which direction?  Whichever one moves in the right y direction.
        direction = -cmp(start.y,end.y)
        # It really can't be zero
        if direction==0:
            raise Exception("??? Trying to draw horizontal line!?")
        (x,y)=((start.x+slope*direction)%self.xmodulus,
               (start.y+direction))
        # non-end-point can't be at x=0, y=0, or y=ymax, so
        # safe to do modulus there.
        while end != (x,y) and x>=0 and y>0 and y<self.ymax:
            rv.append(Point(x,y,self))
            x+=slope*direction    # Slope only affects one coordinate
            x%=self.xmodulus
            y+=direction
        if end != (x,y):
            raise Exception("Could not complete line %s--%s"%
                            (str(start),str(end)))
        return rv

    def pathbetween(self,p1,p2):
        "Return list of points along 'path' between p1 and p2 (which are on one line); i.e. the endpoints and any wraparounds that may be happening."
        slope=self.slopebetween(p1,p2)
        if not slope:
            return []
        direction=cmp(p2.y-p1.y,0)
        rv=[p1]
        # Shouldn't need this, but you know... sometimes you pass a point
        # into the wrong knot or something...
        counter=0
        p2intercept=p2.lines(nowrap=True)[(slope+1)/2]
        while p1 != p2:
            p1intercept=p1.lines(nowrap=True)[(slope+1)/2]
            if p1intercept==p2intercept:
                # They're a simple line apart
                p1=p2
            else:
                # sort of have to do this by condition and not by formula
                # Note that the wraparound points, being one *past* the 
                # wraparound do not have the same y-coordinate.
                if direction*slope > 0:
                    # Moving rightward
                    # !! Note that these are semi-illegal Points, with
                    # x-coordinates outside of [0,self.modulus)
                    rv.append(Point(self.xmodulus+1,
                                    p1.y+slope*(self.xmodulus+1-p1.x), self))
                    p1=Point(-1,
                              p1.y+slope*(self.xmodulus-1-p1.x), self)
                else:
                    rv.append(Point(-1,
                                     p1.y+slope*(-1-p1.x), self))
                    p1=Point(self.xmodulus+1,
                             p1.y+slope*(1-p1.x), self)
            rv.append(p1)
            counter+=1
            if counter>300:
                raise Exception("Seemingly endless path.")
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

    def svgout(self,stroke_width=0.3,scale=20,circle_radius=0.3,
               startat=None,coloriter=None,crossings=True):
        try:
            if type(SVGdraw)!=type(__builtins__):
                return None
        except NameError:
            return None

        cols=['#000000', 
              '#800000', '#808000', '#008080', '#000080',
              '#ff2000', '#ffff20', '#20ffff', '#0020ff',
              '#ff0080', '#ff8000', '#8000ff', '#80ff00']
        svg=SVGdraw.svg(width="%dpx"%((self.xmodulus+2)*scale),
                        height="%dpx"%((self.ymax+2)*scale),
                        viewBox=[-1, -1, self.xmodulus+2,
                                  self.ymax+2])
        defs=SVGdraw.defs(id="defs")
        plusmask=SVGdraw.SVGelement("mask",
                                    attributes={"id":"plusmask"})
        minusmask=SVGdraw.SVGelement("mask",
                                     attributes={"id":"minusmask"})
        r=SVGdraw.rect(x=-1,y=-1,width=self.xmodulus+2,height=self.ymax+2,
                       fill='white')
        plusmask.addElement(r)
        minusmask.addElement(r)
        defs.addElement(plusmask)
        defs.addElement(minusmask)
        svg.addElement(defs)
        maingroup=SVGdraw.group(id="main")
        # I've come to expect them this way up...
        maingroup.attributes['transform']='scale(1,-1) translate(0,%d)'% \
            (-self.ymax)
        svg.addElement(maingroup)
        # Positive slopes and negative slopes.
        plus=SVGdraw.group(id="plus",mask="url(#plusmask)")
        minus=SVGdraw.group(id="minus",mask="url(#minusmask)")
        maingroup.addElement(plus)
        maingroup.addElement(minus)
        circgroup=SVGdraw.group(id="circgroup")
        maingroup.addElement(circgroup)
        strands=self.strands(self.pivots[0])
        circuit=None
        if coloriter is None:
            if len(strands)>1:
                # Multistranded; different color thingie.
                def multicoloriter():
                    counter=0
                    lastcircuit=None
                    while True:
                        if circuit != lastcircuit:
                            lastcircuit=circuit
                            counter+=1
                        yield cols[counter%len(cols)]
                coloriter=multicoloriter()
            else:
                def singlecoloriter(): # for singlestranders!
                    colcounter=0
                    colordiv=len(self.pivots)/6
                    while True:
                        yield cols[int(colcounter/colordiv)%len(cols)]
                        colcounter+=1
                coloriter=singlecoloriter()

            
        for circuit in strands:
            for i in range(0,len(circuit)):
                here=circuit[i]
                nxt=circuit[(i+1)%len(circuit)]
                col=coloriter.next()
                path=self.pathbetween(here,nxt)
                pathstring=""
                for j in range(0,len(path),2):
                    pathstring+=" M %d %d L %d %d"% \
                        (path[j].x,path[j].y,
                         path[j+1].x,path[j+1].y)
                pathelt=SVGdraw.path(pathstring,stroke_width=stroke_width,
                                     stroke=col)
                if self.slopebetween(here,nxt)>0:
                    plus.addElement(pathelt)
                else:
                    minus.addElement(pathelt)
        for i in self.pivots:
            c=SVGdraw.circle(cx=i.x, cy=i.y, r=circle_radius,
                             fill='black')
            circgroup.addElement(c)
        circgroup.addElement(SVGdraw.path("M 0 -1 l 0 %d M %d -1 l 0 %d"% \
                                              (self.ymax+2,self.xmodulus,
                                               self.ymax+2),
                                          stroke='black',
                                          stroke_width=0.03))

        if crossings and len(strands)==1:
            # not handling multistrand with crossings yet...
            circuit=strands[0]
            masked=set()
            over=True
            masks=[minusmask,plusmask]
            # This is a little redundant...
            oncircuit=self.oncircuit(circuit)
            for i in range(0,len(circuit)):
                here=circuit[i]
                nxt=circuit[(i+1)%len(circuit)]
                points=self.pointsbetween(here,nxt)
                slope=self.slopebetween(here,nxt)
                for p in points:
                    if str(p) in masked:
                        over=not over
                        continue # don't do it twice, no need.
                    # only count places that are crossed TWICE.
                    if len(filter(lambda x: x==p, oncircuit))<2:
                        continue
                    if over:
                        mask=masks[(1+slope)/2]
                    else:
                        mask=masks[(1-slope)/2]
                    r=SVGdraw.rect(x=p.x-0.5, y=p.y-0.5,
                                   width=1, height=1,
                                   fill="black",
                                   transform="rotate(45,%d,%d)"%(p.x,p.y))
                    mask.addElement(r)
                    # If it's on the edge, duplicate it on the other side
                    # for ease of viewing.
                    if p.x==0:
                        mask.addElement(SVGdraw.rect(x=self.xmodulus-0.5,
                                                     y=p.y-0.5,
                                                     width=1, height=1,
                                                     fille="#111",
                                                     transform=
                                                     "rotate(45,%d,%d)"%
                                                     (self.xmodulus,p.y)))
                    masked.add(str(p))
                    over=not over
        return svg

def out2file(knot, filename, *args, **kwargs):
    f=open(filename,"w")
    d=SVGdraw.drawing()
    d.svg=knot.svgout(*args,**kwargs)
    f.write(d.toXml())
    f.close
