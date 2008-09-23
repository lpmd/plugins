#
#
#
class LPMDFormat:

  def __init__(self, dc, args={}): pass

  def Show(self): print "*** lpmd plugin information ***"

  def ReadHeader(self, f):
      h = f.readline()
      if not h.startswith('LPMD '): raise self.dc.WrongFormat('lpmd', 'File is not in LPMD format')

  def ReadConfig(self, f):
      config = {'NATOMS': 0, 'AX': 1.0, 'AY': 0.0, 'AZ': 0.0, 'BX': 0.0, 'BY': 1.0, 'BZ': 0.0, 'CX': 0.0, 'CY': 0.0, 'CZ': 1.0}
      line = f.readline()
      if line == '': return None
      config['NATOMS'] = int(line)
      cell = [float(x) for x in f.readline().strip().split()]
      k = 0
      for tag in ['AX', 'AY', 'AZ', 'BX', 'BY', 'BZ', 'CX', 'CY', 'CZ']:
          config[tag] = cell[k]
          k = k + 1
      nats = 0
      while True:
            line = f.readline()
            if line == '': break
            if line.startswith(': END'): break
            if line.startswith(': '):
               lspl = line.split()
               tag, value = lspl[1], lspl[3]
               config[tag] = self.dc.condreducer(value)
            else: 
               lspl = line.split()
               config['CHEMSYMBOL_%d' % nats] = lspl.pop(0)
               config['X_%d' % nats] = float(lspl.pop(0))
               config['Y_%d' % nats] = float(lspl.pop(0))
               config['Z_%d' % nats] = float(lspl.pop(0))
               if len(lspl) > 0: 
                  config['VX_%d' % nats] = float(lspl.pop(0))
                  config['VY_%d' % nats] = float(lspl.pop(0))
                  config['VZ_%d' % nats] = float(lspl.pop(0))
               if len(lspl) > 0: 
                  config['AX_%d' % nats] = float(lspl.pop(0))
                  config['AY_%d' % nats] = float(lspl.pop(0))
                  config['AZ_%d' % nats] = float(lspl.pop(0))
               nats = nats + 1
      return config

  def ReadMany(self, f):
      data = []
      while True:
        c = self.ReadConfig(f)
        if c == None: break
        data.append(c)
      return data

  def WriteHeader(self, f):
      f.write('LPMD 2.0\n')

  def WriteConfig(self, f, config):
      f.write('%d\n' % config['NATOMS'])
      if config.has_key('AX'): 
         f.write('%8.8f  %8.8f  %8.8f  ' % (config['AX'], config['AY'], config['AZ']))
         f.write('%8.8f  %8.8f  %8.8f  ' % (config['BX'], config['BY'], config['BZ']))
         f.write('%8.8f  %8.8f  %8.8f\n' % (config['CX'], config['CY'], config['CZ']))
      else: f.write('1.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0\n')
      tags_used = ['NATOMS', 'AX', 'AY', 'AZ', 'BX', 'BY', 'BZ', 'CX', 'CY', 'CZ']
      for tag in config.keys():
          if not '_' in tag and tag not in tags_used: f.write(': %s = %s\n' % (tag, str(config[tag])))
      for i in range(0, config['NATOMS']):
          tags_used = []
          if config.has_key('CHEMSYMBOL_%d' % i):
             f.write('%s  %8.8f  %8.8f  %8.8f  ' % (config['CHEMSYMBOL_%d' % i], config['X_%d' % i], config['Y_%d' % i], config['Z_%d' % i]))
             tags_used.extend(['CHEMSYMBOL_%d' % i, 'X_%d' % i, 'Y_%d' % i, 'Z_%d' % i])
          if config.has_key('VX_%d' % i):
             f.write('%8.8f  %8.8f  %8.8f  ' % (config['VX_%d' % i], config['VY_%d' % i], config['VZ_%d' % i]))
             tags_used.extend(['VX_%d' % i, 'VY_%d' % i, 'VZ_%d' % i])
          if config.has_key('AX_%d' % i):
             f.write('%8.8f  %8.8f  %8.8f  ' % (config['AX_%d' % i], config['AY_%d' % i], config['AZ_%d' % i]))
             tags_used.extend(['AX_%d' % i, 'AY_%d' % i, 'AZ_%d' % i])
          for tag in config.keys():
              if tag.endswith('_%d' % i) and tag not in tags_used: f.write(': %s = %s\n' % (tag, str(config[tag])))
      f.write(': END\n')

  def WriteMany(self, f, configs):
      for c in configs: self.WriteConfig(f, c)

#
#
#
def pluginLoader(dc, args={}): return LPMDFormat(dc, args)

