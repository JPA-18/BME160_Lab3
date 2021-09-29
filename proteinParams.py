#!/usr/bin/env python3
# Name: Jonathan Postel (jpostel)
# Group Members: Amr Makhamreh (amakhamr)

class ProteinParam :
    ''' This program takes an amino acid sequence of arbitrary length and returns some characterictics and features of the polypeptide
    Input: A string of arbitrary length that should be an amino acid sequence. Other characters are ignored
    Output: The number of amino acids in the string, the molecular weight of the peptide, the molar extinction coefficient, the mass extinctioni coefficient, and the isoelectric point '''
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        ''' This method initializes the objects protein, and aaCompDict when they are instantiated in the class ProteinParam'''
        self.protein = protein.upper()
        aaCompDict ={aa:self.protein.count(aa) for aa in ProteinParam.aa2mw.keys()}  #taking input aa sequence and looking at EACH element and comparing it to the keys in the dictionary provided and if match they get stored
        self.aaCompDict = aaCompDict
    def aaCount (self):
        ''' Counts the number of amino acids in the dictionary that was built from the user's input '''
        return sum(self.aaCompDict.values())


    def pI (self):
        '''calculating specific charge so the pH will be neutral. Needs to be prescise to the 100's
        can limit range to between 0-14
        first test can be 7 and if its positive dont have to check below 7 so new range is 7-14
        next test would be the middle point of the new range and update the range depending if positive or negative and converge to where
        we know when we are done when upper - lower is less than 0.01
        More formally this method finds the isoelectric point of the inputed peptide by calling the _charge_ method. It does so using a binary search instead of checking through every appropriate pH value
        '''
        if (self.aaCount() > 0):
            low = 0.00
            high = 14.00
            while ((high - low) > 0.01):
                cmid = (high + low) / 2
                midInitial = self._charge_(cmid)
                if (midInitial > 0):
                    low = cmid
                else:
                    high = cmid
            else:
                return (cmid)
        else:
            return (0)

    def aaComposition (self):
        return self.aaCompDict

    def _charge_ (self, pH):
        ''' Calculates the charge of the user's input '''
        self.pH = pH
        posCharge = 0
        negCharge = 0
        cTerminus = 10 ** (self.pH) / (10 ** (self.pH) + 10 ** ProteinParam.aaCterm) #calculating charge of c-terminus
        nTerminus = 10 ** (ProteinParam.aaNterm) / ( 10 ** (ProteinParam.aaNterm) + 10 ** (self.pH)) #calculating charge of n-terminus

        #calculates the positive charge
        for aa in self.aaCompDict and ProteinParam.aa2chargePos:
            posCharge += self.aaCompDict[aa] * ( (10 ** (ProteinParam.aa2chargePos[aa])) / (10 ** (ProteinParam.aa2chargePos[aa]) + (10 ** (self.pH))))
        #calculates the negative charge
        for aa in self.aaCompDict and ProteinParam.aa2chargeNeg:
            negCharge += self.aaCompDict[aa] * ((10 ** (self.pH)) / (10 ** (ProteinParam.aa2chargeNeg[aa]) + (10 ** (self.pH))))
        return ((posCharge + nTerminus) - (negCharge + cTerminus))
    def molarExtinction (self):
        ''' Calculates the molar extinction coefficient'''
        molExtinction = 0
        for aa in self.aaCompDict and ProteinParam.aa2abs280:
            molExtinction += self.aaCompDict[aa] * ProteinParam.aa2abs280[aa]
        return molExtinction

    def massExtinction (self):
        ''' Calculates the mass extinction coefficient '''
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight (self):
        ''' calculates the molecular weight of the user's input '''
        aaSum = 0
        for aa in self.protein and ProteinParam.aa2mw:
            if aa in self.aaCompDict.keys():
                aaSum += self.aaCompDict[aa] * (ProteinParam.aa2mw[aa] - ProteinParam.mwH2O)
            else:
                pass

        return ProteinParam.mwH2O + aaSum

import sys
def main():
    ''' write method docstring'''
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")
        myAAcomposition = myParamMaker.aaComposition()
        keys = list(myAAcomposition.keys())
        keys.sort()
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present
        for key in keys :
            print ("\t{} = {:.2%}".format(key, myAAcomposition[key]/myAAnumber))

        inString = input('protein sequence?')


if __name__ == "__main__":
    main()
