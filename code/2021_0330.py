from __future__ import division
import numpy as np
from matplotlib.pyplot import cm
import string
from neuron import h
import numbers
from neuron import gui
from neuron.units import ms, mV
import numpy as np
import matplotlib.pyplot as plt

# a helper library, included with NEURON
h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')
h.celsius = 25


class Cell:
    def __init__(self, name=None, soma=None, apic=None, dend=None, axon=None):
        self.soma = soma if soma is not None else []
        self.apic = apic if apic is not None else []
        self.dend = dend if dend is not None else []
        self.axon = axon if axon is not None else []
        if name is not None:
            self.name = name
        else:
            name = "Retina Ganglion Cell"
        self.all = self.soma + self.apic + self.dend + self.axon

    def delete(self):
        self.soma = None
        self.apic = None
        self.dend = None
        self.axon = None
        self.all = None

    def __str__(self):
        return self.name


def load(filename, fileformat=None, cell=None, use_axon=True, xshift=0, yshift=0, zshift=0):
    """
    Load an SWC from filename and instantiate inside cell. Code kindly provided
    by @ramcdougal.
    Args:
        filename = .swc file containing morphology
        cell = Cell() object. (Default: None, creates new object)
        filename = the filename of the SWC file
        use_axon = include the axon? Default: True (yes)
        xshift, yshift, zshift = use to position the cell
    Returns:
        Cell() object with populated soma, axon, dend, & apic fields
    Minimal example:
        # pull the morphology for the demo from NeuroMorpho.Org
        from PyNeuronToolbox import neuromorphoorg, load
        with open('c91662.swc', 'w') as f:
            f.write(neuromorphoorg.morphology('c91662'))
        cell = load(filename)
    """

    if cell is None:
        # Python 3 不适用，直接改掉了 ————Yuang Cheng
        # cell = Cell(name=string.join(filename.split('.')[:-1]))
        cell = Cell(name="Retina Ganglion Cell")

    if fileformat is None:
        fileformat = filename.split('.')[-1]

    name_form = {1: 'soma[%d]', 2: 'axon[%d]', 3: 'dend[%d]', 4: 'apic[%d]'}

    # load the data. Use Import3d_SWC_read for swc, Import3d_Neurolucida3 for
    # Neurolucida V3, Import3d_MorphML for MorphML (level 1 of NeuroML), or
    # Import3d_Eutectic_read for Eutectic.
    if fileformat == 'swc':
        morph = h.Import3d_SWC_read()
    elif fileformat == 'asc':
        morph = h.Import3d_Neurolucida3()
    else:
        raise Exception('file format `%s` not recognized' % (fileformat))
    morph.input(filename)

    # easiest to instantiate by passing the loaded morphology to the Import3d_GUI
    # tool; with a second argument of 0, it won't display the GUI, but it will allow
    # use of the GUI's features
    i3d = h.Import3d_GUI(morph, 0)

    # get a list of the swc section objects
    swc_secs = i3d.swc.sections
    swc_secs = [swc_secs.object(i) for i in range(int(swc_secs.count()))]

    # initialize the lists of sections
    sec_list = {1: cell.soma, 2: cell.axon, 3: cell.dend, 4: cell.apic}

    # name and create the sections
    real_secs = {}
    for swc_sec in swc_secs:
        cell_part = int(swc_sec.type)

        # skip everything else if it's an axon and we're not supposed to
        # use it... or if is_subsidiary
        if (not (use_axon) and cell_part == 2) or swc_sec.is_subsidiary:
            continue

        # figure out the name of the new section
        if cell_part not in name_form:
            raise Exception('unsupported point type')
        name = name_form[cell_part] % len(sec_list[cell_part])

        # create the section
        sec = h.Section(name=name, cell=cell)

        # connect to parent, if any
        if swc_sec.parentsec is not None:
            sec.connect(real_secs[swc_sec.parentsec.hname()](swc_sec.parentx))

        # define shape
        if swc_sec.first == 1:
            h.pt3dstyle(1, swc_sec.raw.getval(0, 0), swc_sec.raw.getval(1, 0),
                        swc_sec.raw.getval(2, 0), sec=sec)

        j = swc_sec.first
        xx, yy, zz = [swc_sec.raw.getrow(i).c(j) for i in range(3)]
        dd = swc_sec.d.c(j)
        if swc_sec.iscontour_:
            # never happens in SWC files, but can happen in other formats supported
            # by NEURON's Import3D GUI
            raise Exception('Unsupported section style: contour')

        if dd.size() == 1:
            # single point soma; treat as sphere
            x, y, z, d = [dim.x[0] for dim in [xx, yy, zz, dd]]
            for xprime in [x - d / 2., x, x + d / 2.]:
                h.pt3dadd(xprime + xshift, y + yshift, z + zshift, d, sec=sec)
        else:
            for x, y, z, d in zip(xx, yy, zz, dd):
                h.pt3dadd(x + xshift, y + yshift, z + zshift, d, sec=sec)

        # store the section in the appropriate list in the cell and lookup table
        sec_list[cell_part].append(sec)
        real_secs[swc_sec.hname()] = sec

    cell.all = cell.soma + cell.apic + cell.dend + cell.axon
    return cell


Ganglion = load("LY25-RGC10.CNG.swc", 'swc')

Ganglion.soma[0].diam = 15
Ganglion.soma[0].L = 15
Ganglion.soma[1].diam = 15
Ganglion.soma[1].L = 15
Ganglion.soma[2].diam = 15
Ganglion.soma[2].L = 15
Ganglion.soma[0].cm = 6
Ganglion.soma[1].cm = 6
Ganglion.soma[2].cm = 6

for sec in Ganglion.all:
    sec.insert('nargc')
    sec.insert('calrgc')
    sec.insert('canrgc')
    sec.insert('pas')
    sec.insert('kdrgc')
    sec.insert('kargc')
    for seg in sec:
        seg.nargc.gbar = 160
        seg.calrgc.gbar = 0.15
        seg.canrgc.gbar = 0.2
        seg.pas.g = 0.000025
        seg.kargc.gbar = 75
        seg.kdrgc.gbar = 40
        seg.k_ion.ek = -85
        seg.pas.e = -60
        seg.ca_ion.eca = 115
        seg.na_ion.ena = 50

IC = h.IClamp(Ganglion.soma[0](0.5))
IC.delay = 100  # ms
IC.dur = 400  # ms

soma_v = h.Vector()
soma_v.record(Ganglion.soma[0](0.5)._ref_v)
t = h.Vector()
t.record(h._ref_t)
stim_current = h.Vector()
stim_current.record(IC._ref_i)

# 作图
# threshold current
IC.amp = 0.02  # nA
h.v_init = -63 * mV  # initialization voltage
h.dt = 0.025 * ms  # timestep in milliseconds
h.tstop = 800 * ms  # 1e3 millisecond runtime
h.run()  # run simulation for h.tstop
plt.figure()
plt.plot(t.as_numpy(), soma_v.as_numpy(), color='black',
         linewidth='0.5', label="amplitude = 20 pA")
plt.legend(loc='upper right')
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')
plt.xlim((0, h.tstop))
plt.show()

# midrange current
IC.amp = 0.07  # 70pA
h.v_init = -63 * mV  # initialization voltage
h.dt = 0.025 * ms  # timestep in milliseconds
h.tstop = 800 * ms  # 1e3 millisecond runtime
h.run()  # run simulation for h.tstop
plt.figure()
plt.plot(t.as_numpy(), soma_v.as_numpy(), color='black',
         linewidth='0.5', label="amplitude = 70 pA")
plt.legend(loc='upper right')
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')
plt.xlim((0, h.tstop))
# plt.title('amplitude = 70 pA, duration = 400 ms')
plt.show()

# maximum current
IC.amp = 0.150  # 150pA
h.v_init = -63 * mV  # initialization voltage
h.dt = 0.025 * ms  # timestep in milliseconds
h.tstop = 800 * ms  # 1e3 millisecond runtime
h.run()  # run simulation for h.tstop
plt.figure()
plt.plot(t.as_numpy(), soma_v.as_numpy(), color='black',
         linewidth='0.5', label="amplitude = 150 pA")
plt.legend(loc='upper right')
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')
plt.xlim((0, h.tstop))
# plt.title('amplitude = 150 pA')
plt.show()

# Responses to hyperpolarizing currents 超极化电流
# amplitude = -100pA
IC.amp = -0.1
h.v_init = -63 * mV  # initialization voltage
h.dt = 0.025 * ms  # timestep in milliseconds
h.tstop = 800 * ms  # 1e3 millisecond runtime
h.run()  # run simulation for h.tstop
plt.figure()
plt.plot(t.as_numpy(), soma_v.as_numpy(), color='black',
         linewidth='0.5', label="amplitude = -100 pA")
plt.legend(loc='upper right')
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')
plt.ylim(-120, -50)
plt.xlim((0, h.tstop))
# plt.title('amplitude = -100 pA, duration = 400 ms')
plt.show()

# amplitude = -60pA
IC.amp = -0.06
h.v_init = -63 * mV  # initialization voltage
h.dt = 0.025 * ms  # timestep in milliseconds
h.tstop = 800 * ms  # 1e3 millisecond runtime
h.run()  # run simulation for h.tstop
plt.figure()
plt.plot(t.as_numpy(), soma_v.as_numpy(), color='black',
         linewidth='0.5', label="amplitude = -60 pA")
plt.legend(loc='upper right')
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')
plt.ylim(-120, -50)
plt.xlim((0, h.tstop))
# plt.title('amplitude = -60 pA, duration = 400 ms')
plt.show()

# amplitude = -10pA
IC.amp = -0.01
h.v_init = -63 * mV  # initialization voltage
h.dt = 0.025 * ms  # timestep in milliseconds
h.tstop = 800 * ms  # 1e3 millisecond runtime
h.run()  # run simulation for h.tstop
plt.figure()
plt.plot(t.as_numpy(), soma_v.as_numpy(), color='black',
         linewidth='0.5', label="amplitude = -10 pA")
plt.legend(loc='upper right')
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mV)')
plt.ylim(-120, -50)
plt.xlim((0, h.tstop))
# plt.title('amplitude = -10 pA, duration = 400 ms')
plt.show()
