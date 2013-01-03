#!/usr/bin/env python

try:
    import SVGdraw
except ImportError:
    pass

import sys
import math
import itertools

"""
Working solely in the 'diagonal' system from here on!
"""
class Point:
    """A Point in "diagonal" space, i.e. (x,y) is a lattice-point on another
grid which is tilted 45 degrees from this one.  Every point is associated with
some Knot."""
    def __init__(self,x,y,knot):
        if not isinstance(x,int) and isinstance(y,int) and not (x+y)%2:
            raise Exception("Point must be ints and add to an even number")
        self.x=x
        self.y=y
        self.knot=knot

    def __eq__(self,other):
        "Points can equal 2-element lists or tuples as well as other knots."
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
        # Thing is, when adding things to a set, the set machinery requires a
        # hash function to see if it's already there.
        return hash(self.knot) + 1000*self.x + 10*self.y + 10

    @classmethod
    def Pointlist(cls, knot, lst):
        """Take a list of (x,y) tuples and return a list of Points."""
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
    """A Turks' Head Knot, represented by its list of "pivot-points," which are the
places you would put a pin into your mandrell if you were working it.  It's where
the bights turn.

The coordinates used are as if you took the tilted graph paper you normally work
on and drew a grid passing through every lattice point.  This gives a grid twice as
fine-grained as the original one, and we only use those points whose coordinates
add up to an even number."""
    def __init__(self, ptlist):
        if not ptlist:
            return
        ptlist=Point.Pointlist(self,ptlist)
        self.pivots=ptlist[:]   # do I need to make a copy, you think?
        self.valid=True
        self.reconfigure()

    def reconfigure(self):
        "Do all the initial checks and calculations on a knot."
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
        # (b) xmodulus should == len(self.pivots) (True?)
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
by varying choice of pivots within section?
"""
        from fractions import gcd
        total=sum([e[0] for e in layers])
        if any([gcd(total,e[0])<2 for e in layers]):
            raise Exception("Layers must not be prime to total.")
        pivots=[]
        layertable=[]
        heights=[h[1] for h in layers]
        # Will need an iterator or function that can generate the *combinations*
        # in succession. -- itertools.combinations
        # I don't need to consider shifting because "all combinations" includes all
        # shifts of all combinations.
        def assembleknot(combos):
            """assembleknot(combos)
combos is a list [[r00, r01, ... r0k], [r10, r11, ... r1m], ...] where each list of
r's is a combination of elements (numbers): if s_k is gcd(total,n_k), then we are 
choosing n_k/s_k numbers out of range(0,total/s_k).
Return the knot which is the result of repeating those combinations across their 
respective rows.
"""
            # First, the bottom layer.
            l=[(x,0) for x in range(0,2*total,2)]
            for i in range(0,len(combos)):
                sections=gcd(total, layers[i][0])
                for j in range(0,sections):
                    shift=layers[i][1]%2
                    l.extend([(x+2*j*total/sections+shift,layers[i][1]) for x in combos[i]])
            return Knot(l)
        iters=[]
        for layer in layers:
            [number,height]=layer
            sections=gcd(total,number)
            size=total/sections
            howmany=number/sections
            iters.append(itertools.combinations(range(0,2*size,2), howmany))
        bigiter=itertools.product(*iters)
        results=set([])
        for comb in bigiter:
            k=assembleknot(comb)
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
        # Bug or limitation: Something like a 9x3 TH, in which the strands
        # go from one pivot to another and then directly back again, the
        # two points are on BOTH lines.  We only see one, so we wind up
        # just drawing parallel lines instead of crossing ones
        
        # So how about in such a case it depends where you're coming from?
        l1=p1.lines()
        l2=p2.lines()
        if l1==l2:              # If they're BOTH equal, i.e. on both lines
            return cmp(p1.y, p2.y) # can't be zero.
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
               startat=None,coloriter=None,crossings=True,circradius=None,circscale=1):
        # if circradius is some positive number, try to draw a circular(!) diagram
        # circscale is how much to scale the y-dimension by (how thick a circle)
#        try:
#            if type(SVGdraw)!=type(__builtins__):
#		raise Exception("SVGdraw not a module?")
#                return None
#        except NameError:
#	    raise Exception("No SVGDraw found")
#            return None

        cols=['#000000', 
              '#800000', '#808000', '#008080', '#000080',
              '#ff2000', '#ffff20', '#20ffff', '#0020ff',
              '#ff0080', '#ff8000', '#8000ff', '#80ff00']
        if circradius:
            sz=(2*self.ymax*circscale+2+2*circradius)
            svg=SVGdraw.svg(width="%dpx"%(sz*scale), height="%dpx"%(sz*scale),
                            viewBox=[-sz+self.xmodulus/2.0, -sz, 2*sz, 2*sz])
            def transform(x,y):
                # Have to flip it over...
                r=self.ymax*circscale+circradius-y*circscale
                theta=2*math.pi*x/self.xmodulus-math.pi
                return [sz/2+r*math.cos(theta), sz/2+r*math.sin(theta)]
        else:
            svg=SVGdraw.svg(width="%dpx"%((self.xmodulus+2)*scale),
                            height="%dpx"%((self.ymax+2)*scale),
                            viewBox=[-1, -1, self.xmodulus+2,
                                      self.ymax+2])
            def transform(x,y):
                return [x,y]
                        
        defs=SVGdraw.defs(id="defs")
        plusmask=SVGdraw.SVGelement("mask",
                                    attributes={"id":"plusmask"})
        minusmask=SVGdraw.SVGelement("mask",
                                     attributes={"id":"minusmask"})
        if circradius:
            sz=1+2*self.ymax*circscale+2*circradius # Whatever, something big.
            r=SVGdraw.rect(x=-sz, y=-sz, width=sz*2,height=sz*2,fill='white')
        else:
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
                # Multistranded; color it by strand.
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
            # If there's a startat parameter, and it appears in this list,
            # slosh the list around so it's first
            if startat and startat in circuit:
                ind=circuit.index(startat)
                circuit=circuit[ind:]+circuit[0:ind]
            for i in range(0,len(circuit)):
                here=circuit[i]
                nxt=circuit[(i+1)%len(circuit)]
                col=coloriter.next()
                if type(col)==int: # let iterator generate indexes
                    col=cols[col%len(cols)]
                if circradius:
                    path=[here,nxt]
                else:
                    path=self.pathbetween(here,nxt)
                pathstring=""
                for j in range(0,len(path),2):
                    # Had hoped that transform() would have been enough, but we need
                    # to go through all the intermediate lattice-points when doing
                    # circular plots, to curve around in the right direction.
                    if circradius:
                        betweens=self.pointsbetween(path[j],path[j+1])
                        pathstring+=" M %f %f "%tuple(transform(path[j].x,path[j].y))
                        for k in range(0,len(betweens)):
                            pathstring+=" L %f %f "% \
                                tuple(transform(betweens[k].x,betweens[k].y))
                        pathstring+="L %f %f "% \
                            tuple(transform(path[j+1].x, path[j+1].y))
                    else:
                        pathstring+=" M %f %f L %f %f"% \
                            (tuple(transform(path[j].x,path[j].y)+
                                   transform(path[j+1].x,path[j+1].y)))
                pathelt=SVGdraw.path(pathstring,stroke_width=stroke_width,
                                     stroke=col,fill="none")
                if self.slopebetween(here,nxt)>0:
                    plus.addElement(pathelt)
                else:
                    minus.addElement(pathelt)
        for i in self.pivots:
            cr=transform(i.x, i.y)
            c=SVGdraw.circle(cx=cr[0], cy=cr[1], r=circle_radius,
                             fill='black')
            circgroup.addElement(c)
        if not circradius:
            # Mark the wraparound point.
            circgroup.addElement(SVGdraw.path("M 0 -1 l 0 %d M %d -1 l 0 %d"% \
                                                  (self.ymax+2,self.xmodulus,
                                                   self.ymax+2),
                                              stroke='black',
                                              stroke_width=0.03))
        # Somehow I want to *note* when a knot is single-strand or
        # multistrand.
        circgroup.addElement(SVGdraw.text(x=0.2,y=0,
                                          text=str(len(strands)),
                                          fill='#000408',
                                          font_size=1,
                                          font_family='sans-serif',
                                          transform='scale(1,-1)'))

        if crossings:
            # Try multistrand crossings?  (not working right)
            # Need *ALL* the crossing points though.
            oncircuit=[]
            for circuit in strands:
                oncircuit.extend(self.oncircuit(circuit))
            masked=set()
            over=0
            masks=[minusmask,plusmask]
            # How about this?  For each horizontal line _that has intersections on it_,
            # all crossings go in one direction, and that direction alternates.
            #
            # How do we find those lines?
            points=[]
            for circuit in strands:
                for i in range(0,len(circuit)):
                    here=circuit[i]
                    nxt=circuit[(i+1)%len(circuit)]
                    points+=self.pointsbetween(here,nxt)
            heights=[]
            howmanyhits=dict()
            for p in points:
                howmanyhits[p]=howmanyhits.get(p,0)+1
            howmanyhits=[(p,howmanyhits[p]) for p in howmanyhits.keys()]
            howmanyhits=filter((lambda x: x[1]>1), howmanyhits)
            heights=[x[0].y for x in howmanyhits]
            heights.sort()
            # No "sort unique" so just keep track of the last one we saw and skip it.
            # DOESN'T WORK EITHER BUT BETTER THAN BEFORE XXXXXX
            # (testing with python ./knots.py -l 18 17 6 32 6 37)  Works with more
            # symmetrical designs.
            last=None
            for h in heights:
                if h==last:
                    continue
                last=h
                mask=masks[over]
                over=1-over
                for x in range(0,self.xmodulus):
                    p=Point(x,h,self)
                    if p in self.pivots:
                        continue # Skip pivot-points.
                    tp1=transform(p.x-0.5, p.y-0.5)
                    tp2=transform(p.x-0.5, p.y+0.5)
                    tp3=transform(p.x+0.5, p.y+0.5)
                    tp4=transform(p.x+0.5, p.y-0.5)
                    tp=transform(p.x, p.y)
                    if circradius:
                        r=SVGdraw.circle(fill="black",
                                         cx=tp[0], cy=tp[1], r=0.6)
                    else:
                        angle=45 
                        r=SVGdraw.polygon(fill="black",
                                          points=[tp1,tp2,tp3,tp4],
                                          transform="rotate(%f,%f,%f)"% \
                                              (angle, tp[0], tp[1]))
                    mask.addElement(r)
                    # maingroup.addElement(r)
                    # If it's on the edge, duplicate it on the other side
                    # for ease of viewing.
                    if p.x==0 and not circradius:
                        mask.addElement(SVGdraw.rect(x=self.xmodulus-0.5,
                                                     y=p.y-0.5,
                                                     width=1, height=1,
                                                     fill="#111",
                                                     transform=
                                                     "rotate(45,%d,%d)"%
                                                     (self.xmodulus,p.y)))
        return svg

def out2file(knot, filename, *args, **kwargs):
    f=open(filename,"w")
    d=SVGdraw.drawing()
    d.svg=knot.svgout(*args,**kwargs)
    f.write(d.toXml())
    f.close

def usage():
    # Word this better; option args and argv are conflated.
    print """Usage: %s [opts] [args]
\t-h/--help
\t[-n] [-c cols] [-r rad] [-k] -t leads bights
\t[-n] [-c cols] [-r rad] [-k] [-s] [-a] -l n1 h1 n2 h2...
\t[-n] [-c cols] [-r rad] [-k] '[(x1,y1),(x2,y2)...]'

  -h --help:\t\t\tPrint this usage information
  -n --nocrossing:\t\tDon't show crossovers
  -t --turks-head l b:\t\tSimple l x b Turks Head
  -c --colors=num:\t\tNumber of colors to use
  -l --layers n1 h1 n2 h2 ...:\tSearch for "layered" knot
  -s --single:\t\t\tWhen using -l/--layers, only find single-strand knots
  -r --radius=rad:\t\tCircular plot with given inner radius
  -k --knot-only:\t\tJust print out knot (default: output SVG)
  -a --all:\t\t\tUsed with -l; show all knots found; implies -k.
"""%sys.argv[0]


def errorsvg(msg):
    s=SVGdraw.svg(width='5cm', height='1cm')
    s.addElement(SVGdraw.text(0,0,msg))
    d=SVGdraw.drawing()
    d.svg=s
    print d.toXml()
    return s

if __name__=='__main__':
    from getopt import getopt
    (options, argv)=getopt(sys.argv[1:],"hntlc:sr:ka",
                           ["help","nocrossing","turks-head","layers","colors=",
                            "single", "radius=","knot-only","all"])
    opts={e[0]:e[1] for e in options}
    if opts.has_key("-h") or opts.has_key("--help"):
        usage()
        exit(0)
    if opts.has_key("-l") or opts.has_key("--layers"):
        l=[]
        if len(argv)%2:
            print "Layers requires an even number of arguments (num, height...)"
            exit(1)
        for i in range(0,len(argv)-1,2):
            l.append((int(argv[i]),int(argv[i+1])))
        possibles=Knot.Layers(l)
        if not possibles:
            # Which is better?
            # print "None found"
            errorsvg("None found")
            exit(1)
        elif opts.has_key("-s") or opts.has_key("--single"):
            if all([len(k.strands())>1 for k in possibles]):
                print "Only multistrand knots found."
                exit(2)
        if opts.has_key("-a") or opts.has_key("--all"):
            print str(possibles)
            exit(0)
        k=possibles.pop()
    elif opts.has_key('-t') or opts.has_key('--turks-head'):
        try:
            k=Knot.TH(int(argv[0]),int(argv[1]))
        except:
            print "Simple turks-head, specify <leads> <bights>"
            exit(1)
    else:
        import json
        try:
            s="".join(argv).replace('(','[').replace(')',']')
            l=json.loads(s)
        except ValueError:
            usage()
            exit(1)
        k=Knot(l)
    citer=None
    if opts.has_key('-c') or opts.has_key("--colors"):
        rep=float(opts.get('-c') or opts['--colors'])
        per=len(k.pivots)/rep
        if rep>0:
            def citergen():
                count=0
                while True:
                    yield int(count/per)
                    count+=1
            citer=citergen()
    if opts.has_key("-r") or opts.has_key("--radius"):
        rad=float(opts.get("-r") or opts["--radius"])
    else:
        rad=None
    if opts.has_key("-k") or opts.has_key("--knot-only"):
        # Don't output the svg, just the knot.
        print str(k)
        exit(0)
    s=k.svgout(crossings=not (opts.has_key('-n') or 
                              opts.has_key('--nocrossings')),
               coloriter=citer,circradius=rad)
    d=SVGdraw.drawing()
    d.svg=s
    print d.toXml()

# An irregular one I found experimenting:
# [(1,1),(2,4),(3,1),(4,16),(5,1),(6,10),(6,4),(7,1),(9,9),(9,1),(10,4),(11,1),(13,9),(13,1),(15,1),(16,8),(17,1),(19,1),(21,11),(21,7),(21,1),(22,4),(23,1),(25,1),(26,4),(27,1),(28,16),(28,12),(29,1),(30,6),(30,4),(31,1)]
# (16 wide, but make sure it's tall enough!)

# An interesting conundrum: 
# [(0,2),(0,0),(4,2),(4,0),(5,17),(5,13),(5,9),(5,5),(8,2),(8,0),(11,17),(11,13),(11,9),(11,5),(12,2),(12,0)]
# Irregular bottom.  Length of bottom != wrap zone (I think).  Says it's 
# not tyable, but is it?  If I set the xmodulus to 16 by hand, it's fine;
# generates an SVG and all, tyable in 4 strands.  But as it is, the program
# sets the xmodulus to 14.
# XXXXX But its svg file looks like it has bad crossings!
 
# [(0,0),(1,9),(2,4),(3,13),(4,0),(5,9),(6,4),(7,13),(8,0),(9,9),(10,4),(11,13),(12,0),(13,9),(14,4),(15,13)] is okay though.

# Patrick's:
# [(0,0),(1,13),(2,0),(3,17),(4,0),(5,13),(6,0),(7,33),(8,0),(9,13),(10,0),(11,17),(12,0),(13,13),(14,0),(15,33),(16,0),(17,13),(18,0),(19,17),(20,0),(21,13),(22,0),(23,33),(24,0),(25,13),(26,0),(27,17),(28,0),(29,13),(30,0),(31,33),(32,0),(33,13),(34,0),(35,17),(36,0),(37,13),(38,0),(39,33),(40,0),(41,13),(42,0),(43,17),(44,0),(45,13),(46,0),(47,33)]
