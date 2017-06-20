from __future__ import print_function
from pfnet import CustomParser, Network, ParserError

class DummyParser(CustomParser):

    def init(self):
        self.parsed_data = {}

    def parse(self,filename,num_periods=1):       
 
        if filename.split('.')[-1] != 'dummy':
            raise ParserError('invalid extension')

        net = Network(num_periods)

        return net

    def set(self,key,value):
        pass

    def show(self):
        print('this is a dummy parser')

    def write(self,net,filename):
        pass
