import os
import unittest
from subprocess import call, STDOUT

class TestExamples(unittest.TestCase):

    def test_constraints(self):
        
        FNULL = open(os.devnull, 'w')
        retcode = call(["python", "./examples/constraints.py", "./data/ieee14.mat"],
                       stdout=FNULL,
                       stderr=STDOUT)
        
        self.assertEqual(retcode, 0)

    def test_contingencies(self):
        
        FNULL = open(os.devnull, 'w')
        retcode = call(["python", "./examples/contingencies.py", "./data/ieee14.mat"],
                       stdout=FNULL,
                       stderr=STDOUT)
        
        self.assertEqual(retcode, 0)

    def test_functions(self):
        
        FNULL = open(os.devnull, 'w')
        retcode = call(["python", "./examples/functions.py", "./data/ieee14.mat"],
                       stdout=FNULL,
                       stderr=STDOUT)
        
        self.assertEqual(retcode, 0)

    def test_muti_period(self):
        
        FNULL = open(os.devnull, 'w')
        retcode = call(["python", "./examples/multi_period.py", "./data/ieee14.mat"],
                       stdout=FNULL,
                       stderr=STDOUT)
        
        self.assertEqual(retcode, 0)

    def test_networks(self):
        
        FNULL = open(os.devnull, 'w')
        retcode = call(["python", "./examples/networks.py", "./data/ieee14.mat"],
                       stdout=FNULL,
                       stderr=STDOUT)
        
        self.assertEqual(retcode, 0)

    def test_parsers(self):
        
        FNULL = open(os.devnull, 'w')
        retcode = call(["python", "./examples/parsers.py", "./data/ieee14.mat"],
                       stdout=FNULL,
                       stderr=STDOUT)
        
        self.assertEqual(retcode, 0)

    def test_problems(self):
        
        FNULL = open(os.devnull, 'w')
        retcode = call(["python", "./examples/problems.py", "./data/ieee14.mat"],
                       stdout=FNULL,
                       stderr=STDOUT)
        
        self.assertEqual(retcode, 0)

    def test_projections(self):
        
        FNULL = open(os.devnull, 'w')
        retcode = call(["python", "./examples/projections.py", "./data/ieee14.mat"],
                       stdout=FNULL,
                       stderr=STDOUT)
        
        self.assertEqual(retcode, 0)

    def test_start(self):
        
        FNULL = open(os.devnull, 'w')
        retcode = call(["python", "./examples/start.py", "./data/ieee14.mat"],
                       stdout=FNULL,
                       stderr=STDOUT)
        
        self.assertEqual(retcode, 0)

    def test_variables(self):
        
        FNULL = open(os.devnull, 'w')
        retcode = call(["python", "./examples/variables.py", "./data/ieee14.mat"],
                       stdout=FNULL,
                       stderr=STDOUT)
        
        self.assertEqual(retcode, 0)
