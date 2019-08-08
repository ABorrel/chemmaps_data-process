from pydpi.drug import constitution, topology, connectivity, kappa, bcut, basak, estate, moran, moreaubroto, geary, \
    charge, molproperty, moe, fingerprint
from pydpi import pydrug
from molvs import standardize_smiles, Standardizer


from rdkit.Chem.SaltRemover import SaltRemover
from rdkit import Chem
from rdkit.Chem import AllChem

from copy import deepcopy
from os import path, getcwd, remove, system, listdir

import toolbox
import pathFolder
import runExternalSoft
import descriptors3D


LSALTDEF="[Cl,Br,I]\n[Li,Na,K,Ca,Mg]\n[O,N]\n[N](=O)(O)O\n[P](=O)(O)(O)O\n[P](F)(F)(F)(F)(F)F\n[S](=O)(=O)(O)O\n[CH3][S](=O)(=O)(O)\nc1cc([CH3])ccc1[S](=O)(=O)(O)\n[CH3]C(=O)O\nFC(F)(F)C(=O)O\nOC(=O)C=CC(=O)O\nOC(=O)C(=O)O\nOC(=O)C(O)C(O)C(=O)O\nC1CCCCC1[NH]C1CCCCC1\n"
LSALT="[Co]\n[Cd+2]\n"


LSMILESREMOVE=["[C-]#N", "[Al+3]", "[Gd+3]", "[Pt+2]", "[Au+3]", "[Bi+3]", "[Al]", "[Si+4]", "[Fe]", "[Zn]", "[Fe+2]",
               "[Ru+8]", "[Fe+]", "[Sr++]", "[Fe+3]", "[O--]", "[OH-]", "[Mn++]", "[La+3]", "[Lu+3]", "[SH-]", "[Pt+4]",
               "[Fe++]", "[W]", "[Cu+2]", "[Cr+3]", "[Tc+7]", "[Xe]", "[Tl+]", "[Zn+2]", "[F-]", "[C]", "[He]", "N#N",
               "O=O", "Cl[Ra]Cl", "[Mn+2]", "N#[N+][O-]", "II", "[Ga+3]", "[Mo+10]", "[Zn]", "[Fe]", "[Si+4]", "[Al]",
               "[B+3]"]


LKAPA = ['kappa1', 'kappa2', 'kappa3', 'kappam1', 'kappam2', 'kappam3', 'phi']
LBUCUT =["bcutp16","bcutp15","bcutp14","bcutp13","bcutp12","bcutp11","bcutp10",
        "bcutp9","bcutp8","bcutp7","bcutp6","bcutp5","bcutp4","bcutp3",
        "bcutp2","bcutp1"]
LESTATE = ['Smax38', 'Smax39', 'Smax34', 'Smax35', 'Smax36', 'S43', 'Smax30', 'Smax31', 'Smax32', 'Smax33', 'S57',
           'S56', 'S55', 'S54', 'S53', 'S52', 'S51', 'S50', 'Smin49', 'S59', 'S58', 'Smin69', 'Smin68', 'Smin27',
           'Sfinger30', 'Sfinger31', 'Sfinger32', 'Sfinger33', 'Sfinger34', 'Sfinger35', 'Sfinger36', 'Sfinger37',
           'Sfinger38', 'Sfinger39', 'Smax2', 'Smax3', 'Smax4', 'Smax5', 'Smax6', 'Smax7', 'Smin77', 'Smax29', 'Smax37',
           'Smax23', 'Smax22', 'Smax21', 'Smax20', 'Smax27', 'Smax26', 'Smax25', 'Smax24', 'S44', 'S45', 'S46', 'S47',
           'S40', 'S41', 'S42', 'S17', 'Smin44', 'S48', 'S49', 'Smin8', 'Smin29', 'Smin28', 'Sfinger45', 'Sfinger44',
           'Sfinger47', 'Sfinger46', 'Sfinger41', 'Sfinger40', 'Sfinger43', 'Sfinger42', 'Smax47', 'Smin73', 'Smin70',
           'Smin71', 'Sfinger49', 'Sfinger48', 'Smin74', 'Smin75', 'Smin67', 'Smin6', 'Smin9', 'Smin7', 'Smin47',
           'Smax41', 'S79', 'S78', 'Smin19', 'Smax58', 'Smax59', 'S71', 'S70', 'S73', 'S72', 'S75', 'S74', 'S77',
           'S76', 'Smax73', 'Smin78', 'Sfinger56', 'Sfinger57', 'Sfinger54', 'Sfinger55', 'Sfinger52', 'Sfinger53',
           'Sfinger50', 'Sfinger51', 'Smin61', 'Smin60', 'Smin63', 'Smin62', 'Smin65', 'Smin64', 'Sfinger58',
           'Sfinger59', 'Smin48', 'Smin42', 'Smin76', 'Smin41', 'Smin72', 'Smax40', 'Smin40', 'Smax49', 'Smax48',
           'S68', 'S69', 'S66', 'S67', 'S64', 'S65', 'S62', 'S63', 'S60', 'S61', 'Smin54', 'Smax52', 'Sfinger69',
           'Sfinger68', 'Smin50', 'Smin51', 'Smin52', 'Smin53', 'Sfinger63', 'Sfinger62', 'Sfinger61', 'Sfinger60',
           'Sfinger67', 'S10', 'Sfinger65', 'Sfinger64', 'S13', 'S12', 'Sfinger76', 'Smin56', 'S9', 'S8', 'S3', 'S2',
           'S1', 'Smin55', 'S7', 'S6', 'S5', 'S4', 'Smax78', 'Smax45', 'Smax11', 'Sfinger72', 'Smin66', 'Smax44',
           'Smax70', 'Smax71', 'Smax72', 'S14', 'Smax74', 'Smax75', 'Smax76', 'Smax77', 'Smin43', 'Smax8', 'S19',
           'S18', 'Sfinger78', 'Sfinger79', 'Smin45', 'Smax9', 'Sfinger74', 'Sfinger75', 'S11', 'Sfinger77',
           'Sfinger70', 'Sfinger71', 'S15', 'Sfinger73', 'Smax43', 'Smin16', 'Smax42', 'Smax53', 'Smax66', 'Smax65',
           'Smax64', 'Smax63', 'Smax62', 'Smax61', 'Smax60', 'Smin26', 'Smax69', 'Smax68', 'Smax0', 'Smin57', 'Smax1',
           'Smin17', 'Smin36', 'Smin37', 'Smin34', 'Smin35', 'Smin32', 'Smin33', 'Smin30', 'Smin31', 'Smax67', 'Smin46',
           'Smax51', 'Smin38', 'Smin39', 'Smax12', 'Smax13', 'Smax10', 'S16', 'Smax16', 'Smax17', 'Smax14', 'Smax15',
           'Smin20', 'Smax18', 'Smax19', 'Sfinger66', 'Smax56', 'Smax28', 'Smax57', 'Smax54', 'Smin58', 'Smax55', 'S39',
           'S38', 'Smax46', 'S35', 'S34', 'S37', 'S36', 'S31', 'S30', 'S33', 'S32', 'Smin25', 'Smin24', 'Sfinger18',
           'Sfinger19', 'Smin21', 'Smax50', 'Smin23', 'Smin22', 'Sfinger12', 'Sfinger13', 'Sfinger10', 'Sfinger11',
           'Sfinger16', 'Sfinger17', 'Sfinger14', 'Sfinger15', 'Sfinger8', 'Sfinger9', 'Smin4', 'Smin5', 'Smin2',
           'Smin3', 'Smin0', 'Smin1', 'Sfinger1', 'Sfinger2', 'Sfinger3', 'Sfinger4', 'Sfinger5', 'Sfinger6',
           'Sfinger7', 'S22', 'S23', 'S20', 'S21', 'S26', 'S27', 'S24', 'S25', 'Smin59', 'S28', 'S29', 'Smin18',
           'Smin10', 'Smin11', 'Smin12', 'Smin13', 'Smin14', 'Smin15', 'Sfinger29', 'Sfinger28', 'Sfinger27',
           'Sfinger26', 'Sfinger25', 'Sfinger24', 'Sfinger23', 'Sfinger22', 'Sfinger21', 'Sfinger20']
LMOREAUBROTO = ['ATSe1', 'ATSe2', 'ATSe3', 'ATSe4', 'ATSe5', 'ATSe6', 'ATSe7', 'ATSe8', 'ATSp8', 'ATSp3', 'ATSv8',
                'ATSp1', 'ATSp7', 'ATSp6', 'ATSp5', 'ATSp4', 'ATSv1', 'ATSp2', 'ATSv3', 'ATSv2', 'ATSv5', 'ATSv4',
                'ATSv7', 'ATSv6', 'ATSm8', 'ATSm1', 'ATSm2', 'ATSm3', 'ATSm4', 'ATSm5', 'ATSm6', 'ATSm7']
LMORAN = ['MATSv8', 'MATSp4', 'MATSp8', 'MATSv1', 'MATSp6', 'MATSv3', 'MATSv2', 'MATSv5', 'MATSv4', 'MATSv7', 'MATSv6',
          'MATSm8', 'MATSp1', 'MATSm4', 'MATSm5', 'MATSm6', 'MATSm7', 'MATSm1', 'MATSm2', 'MATSm3', 'MATSe4', 'MATSe5',
          'MATSe6', 'MATSe7', 'MATSe1', 'MATSe2', 'MATSe3', 'MATSe8', 'MATSp3', 'MATSp7', 'MATSp5', 'MATSp2']
LGEARY = ['GATSp8', 'GATSv3', 'GATSv2', 'GATSv1', 'GATSp6', 'GATSv7', 'GATSv6', 'GATSv5', 'GATSv4', 'GATSe2', 'GATSe3',
          'GATSv8', 'GATSe6', 'GATSe7', 'GATSe4', 'GATSe5', 'GATSp5', 'GATSp4', 'GATSp7', 'GATSe1', 'GATSp1', 'GATSp3',
          'GATSp2', 'GATSe8', 'GATSm2', 'GATSm3', 'GATSm1', 'GATSm6', 'GATSm7', 'GATSm4', 'GATSm5', 'GATSm8']
LMOE = ['EstateVSA8', 'EstateVSA9', 'EstateVSA4', 'EstateVSA5', 'EstateVSA6', 'EstateVSA7', 'EstateVSA0', 'EstateVSA1',
        'EstateVSA2', 'EstateVSA3', 'PEOEVSA13', 'PEOEVSA12', 'PEOEVSA11', 'PEOEVSA10', 'MTPSA', 'VSAEstate0',
        'VSAEstate1', 'VSAEstate2', 'VSAEstate3', 'VSAEstate4', 'VSAEstate5', 'VSAEstate6', 'VSAEstate7', 'VSAEstate8',
        'LabuteASA', 'PEOEVSA3', 'PEOEVSA2', 'PEOEVSA1', 'PEOEVSA0', 'PEOEVSA7', 'PEOEVSA6', 'PEOEVSA5', 'PEOEVSA4',
        'MRVSA5', 'MRVSA4', 'PEOEVSA9', 'PEOEVSA8', 'MRVSA1', 'MRVSA0', 'MRVSA3', 'MRVSA2', 'MRVSA9', 'slogPVSA10',
        'slogPVSA11', 'MRVSA8', 'MRVSA7', 'MRVSA6', 'EstateVSA10', 'slogPVSA2', 'slogPVSA3', 'slogPVSA0', 'slogPVSA1',
        'slogPVSA6', 'slogPVSA7', 'slogPVSA4', 'slogPVSA5', 'slogPVSA8', 'slogPVSA9', 'VSAEstate9', 'VSAEstate10']

l3D = ['RDFC6', 'MoRSEN11', 'RDFU8', 'RDFU9', 'RDFU2', 'RDFU3', 'MoRSEN5', 'RDFU1', 'RDFU6', 'RDFU7', 'RDFU4', 'RDFU5',
       'Harary3D', 'P2u', 'MoRSEM6', 'MoRSEM7', 'MoRSEM4', 'MoRSEM5', 'MoRSEM2', 'MoRSEM3', 'MoRSEE30', 'MoRSEM1',
       'MoRSEN4', 'MoRSEM8', 'MoRSEM9', 'MoRSEU10', 'MoRSEU11', 'MoRSEU12', 'MoRSEU13', 'MoRSEU14', 'MoRSEU15',
       'MoRSEU16', 'MoRSEU17', 'MoRSEU18', 'MoRSEU19', 'FPSA3', 'FPSA2', 'FPSA1', 'GeDi', 'MoRSEV19', 'MoRSEN10',
       'MoRSEV13', 'SPAN', 'MoRSEV11', 'MoRSEV10', 'MoRSEV17', 'MoRSEV16', 'MoRSEV15', 'MoRSEN16', 'RDFM14', 'RDFM15',
       'RDFM16', 'RDFE9', 'RDFM10', 'RDFM11', 'RDFM12', 'RASA', 'RDFE2', 'RDFE3', 'MoRSEC25', 'MoRSEC24', 'RDFM18',
       'RDFM19', 'grav', 'RDFE5', 'WNSA1', 'WNSA2', 'WNSA3', 'L2p', 'RDFP15', 'RDFP14', 'RDFP17', 'RDFP16', 'RDFP11',
       'RDFP10', 'RDFP13', 'RDFP12', 'MoRSEP30', 'RDFP19', 'RDFP18', 'E2p', 'Dm', 'P3e', 'MoRSEM18', 'MoRSEM19',
       'Petitj3D', 'MoRSEM10', 'MoRSEM11', 'MoRSEM12', 'MoRSEM13', 'MoRSEM14', 'MoRSEM15', 'MoRSEM16', 'MoRSEM17',
       'RDFC27', 'RDFC26', 'RDFC25', 'RDFC24', 'RDFC23', 'RDFC22', 'RDFC21', 'RDFC20', 'MoRSEC28', 'RDFU30', 'RDFC29',
       'RDFC28', 'MoRSEU30', 'L1u', 'L1v', 'L2v', 'L1p', 'RDFP5', 'RDFP4', 'RDFP7', 'RDFP6', 'RDFP1', 'RDFP3', 'RDFP2',
       'L1e', 'RDFP9', 'RDFP8', 'MoRSEP5', 'P1e', 'MoRSEP4', 'PSA', 'MoRSEP7', 'P1p', 'MoRSEP6', 'RDFE18', 'RDFE19',
       'P1v', 'RDFE14', 'RDFE15', 'RDFE16', 'RDFE17', 'PPSA1', 'RDFE11', 'RDFE12', 'PPSA2', 'MoRSEP11', 'MoRSEP10',
       'MoRSEP13', 'MoRSEP12', 'RPCS', 'MoRSEP14', 'MoRSEN9', 'MoRSEN8', 'DPSA1', 'MoRSEC30', 'DPSA3', 'DPSA2',
       'MoRSEN3', 'MoRSEN2', 'MoRSEN1', 'RDFP30', 'E2e', 'MoRSEN17', 'L3e', 'TASA', 'RDFC19', 'MoRSEV14', 'MoRSEM30',
       'MoRSEP8', 'L3v', 'RDFC16', 'L3u', 'RDFV30', 'L3p', 'RDFC14', 'W3DH', 'RDFC15', 'MoRSEC23', 'MoRSEN15',
       'MoRSEP16', 'RPSA', 'P3m', 'MEcc', 'MoRSEC22', 'MoRSEN14', 'MoRSEP1', 'MoRSEN23', 'P3p', 'P3v', 'MoRSEP19',
       'P3u', 'RDFV7', 'RDFC18', 'RDFV6', 'FNSA1', 'RDFC17', 'FNSA3', 'FNSA2', 'RDFC12', 'RDFC13', 'RDFC10', 'RDFC11',
       'P2p', 'RDFV4', 'MoRSEP22', 'RDFV3', 'MoRSEP18', 'RDFV2', 'RDFU21', 'RDFU20', 'RDFU23', 'RDFU22', 'RDFU25',
       'RDFM4', 'RDFU27', 'RDFU26', 'RDFU29', 'RDFU28', 'MoRSEN28', 'MoRSEN29', 'RDFV19', 'RDFV18', 'WPSA2', 'RDFV16',
       'RDFV15', 'RDFV14', 'RDFV13', 'RDFV12', 'RDFV11', 'RDFV10', 'MoRSEN26', 'MoRSEP23', 'MoRSEN27', 'MoRSEN24',
       'MoRSEP25', 'MoRSEN25', 'MoRSEP26', 'MoRSEE8', 'MoRSEE9', 'MoRSEE6', 'MoRSEN22', 'MoRSEE4', 'MoRSEP27',
       'MoRSEE2', 'MoRSEE3', 'RDFU24', 'MoRSEN20', 'MoRSEC9', 'ASPAN', 'RDFE10', 'MoRSEN21', 'Te', 'Vm', 'Vp',
       'MoRSEV18', 'PPSA3', 'Vv', 'RDFE13', 'E2u', 'RDFC30', 'E2v', 'P1m', 'MoRSEV12', 'MoRSEP15', 'MoRSEP17',
       'MoRSEU8', 'MoRSEU9', 'MoRSEU6', 'MoRSEU7', 'MoRSEU4', 'MoRSEU5', 'MoRSEU2', 'MoRSEU3', 'MoRSEU1', 'RDFM29',
       'RDFM28', 'MoRSEN6', 'RDFM21', 'RDFM20', 'RDFM23', 'RDFM22', 'RDFM25', 'RDFM24', 'RDFM27', 'RDFM26', 'RDFM2',
       'RDFM3', 'RDFV5', 'RDFM1', 'RDFM6', 'RDFM7', 'RDFV1', 'RDFM5', 'MoRSEP20', 'MoRSEP21', 'RDFM8', 'RDFM9',
       'MoRSEP24', 'RDFE8', 'RDFV9', 'RDFV8', 'MoRSEC8', 'RDFV17', 'RDFM17', 'WPSA3', 'AGDD', 'MoRSEC1', 'MoRSEC2',
       'MoRSEC3', 'MoRSEC4', 'MoRSEC5', 'MoRSEC6', 'MoRSEC7', 'MoRSEE29', 'MoRSEE28', 'WPSA1', 'MoRSEC29', 'MoRSEE21',
       'MoRSEE20', 'MoRSEE23', 'RDFM13', 'MoRSEE25', 'MoRSEE24', 'MoRSEE27', 'MoRSEE26', 'MoRSEC27', 'MoRSEC26', 'Ae',
       'RDFE1', 'RDFE6', 'RDFE7', 'RDFE4', 'MoRSEV28', 'MoRSEV29', 'MoRSEV26', 'MoRSEV27', 'MoRSEV24', 'MoRSEC20',
       'MoRSEV22', 'MoRSEV23', 'MoRSEV20', 'MoRSEV21', 'MoRSEC12', 'MoRSEC13', 'MoRSEC10', 'MoRSEC11', 'MoRSEC16',
       'MoRSEC17', 'MoRSEC14', 'MoRSEC15', 'P2m', 'MoRSEC18', 'MoRSEC19', 'RDFE30', 'RDFE21', 'RDFE20', 'RDFE23',
       'RDFE22', 'RDFE25', 'RDFE24', 'RDFE27', 'RDFE26', 'RDFE29', 'RDFE28', 'RDFP20', 'RDFP21', 'RDFP22', 'RDFP23',
       'RDFP24', 'RDFP25', 'RDFP26', 'RDFP27', 'RDFP28', 'RDFP29', 'MoRSEV6', 'MoRSEN13', 'MoRSEV30', 'Dv', 'RDFV26',
       'RDFV27', 'RDFV24', 'RDFV25', 'RDFV22', 'RDFV23', 'RDFV20', 'RDFV21', 'L2e', 'MoRSEV5', 'RDFV28', 'RDFV29',
       'MoRSEP3', 'P1u', 'rygr', 'Ve', 'MoRSEE7', 'MoRSEV4', 'MoRSEP2', 'FrTATP', 'MoRSEE5', 'P2v', 'ASA', 'MoRSEC21',
       'MoRSEV3', 'MoRSEE1', 'E3v', 'E3u', 'Ke', 'E3p', 'Km', 'MoRSEV7', 'E3e', 'Kp', 'Kv', 'Ku', 'MoRSEV2', 'RDFC8',
       'E3m', 'RDFC9', 'MoRSEM21', 'MoRSEM20', 'MoRSEM23', 'MoRSEM22', 'MoRSEM25', 'MoRSEM24', 'MoRSEM27', 'MoRSEM26',
       'MoRSEM29', 'MoRSEM28', 'MoRSEN19', 'MoRSEN18', 'L2m', 'MoRSEV9', 'MoRSEV8', 'MoRSEU29', 'MoRSEU28', 'L2u',
       'MoRSEV1', 'MoRSEN12', 'RDFC5', 'MoRSEU21', 'MoRSEU20', 'MoRSEU23', 'MoRSEU22', 'MoRSEU25', 'MoRSEU24',
       'MoRSEU27', 'MoRSEU26', 'RDFC7', 'MoRSEV25', 'Tv', 'Am', 'SEig', 'Tu', 'Tp', 'RDFC1', 'Tm', 'RDFC2', 'PNSA3',
       'PNSA2', 'PNSA1', 'RDFC3', 'MSA', 'MoRSEE22', 'L1m', 'De', 'MoRSEE10', 'MoRSEE11', 'MoRSEE12', 'MoRSEE13',
       'MoRSEP9', 'MoRSEE15', 'MoRSEE16', 'MoRSEE17', 'MoRSEE18', 'MoRSEE19', 'Du', 'Dp', 'Vu', 'P2e', 'E1p', 'E1u',
       'E1v', 'E1m', 'RNCS', 'MoRSEP28', 'MoRSEN30', 'E1e', 'MoRSEN7', 'W3D', 'RDFU18', 'RDFU19', 'RDFU14', 'RDFU15',
       'RDFU16', 'RDFU17', 'RDFU10', 'RDFU11', 'RDFU12', 'RDFU13', 'Ap', 'Au', 'RDFC4', 'MoRSEP29', 'Av', 'L3m',
       'RDFM30', 'MoRSEE14', 'E2m']


loader = pydrug.PyDrug()



def normalize(mol, lout):
    s = Standardizer()
    molstandardized = s.standardize(mol)
    #print molstandardized
    lout.append(molstandardized)



def getLdesc (typeDesc):

    lout = []
    if typeDesc == "1D2D":
        # listdesc
        lout = lout + constitution._constitutional.keys() + ["nheavy"] + molproperty.MolecularProperty.keys() + topology._Topology.keys() +connectivity._connectivity.keys() + LKAPA + LBUCUT + basak._basak.keys() + LESTATE + LMOREAUBROTO + LMORAN + LGEARY + charge._Charge.keys() + LMOE

    if typeDesc == "3D":
        lout = lout + l3D

    return lout



class chemical:

    def __init__(self, name, smiles, psdf=""):

        self.name = name
        self.smi = smiles
        self.psdf = psdf
        self.log = "Init => " + str(smiles) + "\n"
        self.err = 0
        self.inchikey = ""
        # generate smi?

        #smile = runExternalSoft.babelConvertSDFtoSMILE(self.compound["sdf"])
        #self.compound["SMILES"] = smile
            # print smile
        #except:


    def prepareChem(self, prSMIclean):


        psmiclean = prSMIclean + self.name + ".smi"
        if self.name == "":
            self.err = 1
            return 1

        # try if existing
        if path.exists(psmiclean):
            #print psmiclean
            fsmiclean = open(psmiclean, "r")
            llsmiclean = fsmiclean.readlines()
            fsmiclean.close()

            #print "Path exists", psmiclean

            try:
                lsmiclean = llsmiclean[0].strip().split("\t")
                smiclean = lsmiclean[0]
            except:
                smiclean = "ERROR"

            if smiclean == "ERROR":
                self.err = 1
                self.log = self.log + "Normalize SMILES: file output ERROR\n"

            else:
                if len(lsmiclean) < 3:
                    smiclean == "ERROR"
                    self.err = 1
                    return 1
                else:
                    self.smiclean = smiclean
                    self.inchikey = lsmiclean[1]
                    self.psmiclean = psmiclean
                    self.log = self.log + "Prep SMI :" + str(self.smi) + "\n"
                    self.log = self.log + "Prepared SMI :" + str(self.smiclean) + "\n"

        else:
            s = Standardizer()
            mol = Chem.MolFromSmiles(self.smi)
            try:
                out = toolbox.timeFunction(normalize, mol)
                if out == "ERROR":
                    self.err = 1
                    self.log = self.log + "Normalize SMILES: ERROR DURING THE PROCESS\n"
                    fsmiclean = open(psmiclean, "w")
                    fsmiclean.write("ERROR")
                    fsmiclean.close()
                else:
                    molstandardized = out
            except:
                self.err = 1
                self.log = self.log + "Normalize SMILES: ERROR INPUT SMI\n"
                fsmiclean = open(psmiclean, "w")
                fsmiclean.write("ERROR")
                fsmiclean.close()


            if "molstandardized" in locals():

                smilestandadized = Chem.MolToSmiles(molstandardized)

                # remove salt
                # 1.default
                remover = SaltRemover(defnFilename="Salts.txt")
                mol = Chem.MolFromSmiles(smilestandadized)
                molcleandefault = remover(mol)
                # 2. Personal remover
                homeremover = SaltRemover(defnData=LSALT)
                molclean = homeremover(molcleandefault)
                smilesclean = Chem.MolToSmiles(molclean)
                # 3. SMILES remove other manual salts + fragments -> for fragment take one if exactly same compound
                lelem = smilesclean.split(".")
                if len(lelem) > 1:
                    # reduce double, case of several salts are included - 255
                    lelem = list(set(lelem))
                    for smilesdel in LSMILESREMOVE:
                        if smilesdel in lelem:
                            lelem.remove(smilesdel)
                    try:
                        lelem.remove("")  # case of bad smile
                    except:
                        pass
                    if len(lelem) == 1:
                        smilesclean = str(lelem[0])
                    else:
                        # 4. Fragments
                        # Case of fragment -> stock in log file, check after to control
                        self.log = self.log + "Fragments after standardization: " + smilesclean + "\n"
                        self.err = 1
                        fsmiclean = open(psmiclean, "w")
                        fsmiclean.write("ERROR")
                        fsmiclean.close()
                        return 1

                if smilesclean == "":
                    self.log = self.log + "SMILES empty after preparation\n"
                    self.err = 1
                    fsmiclean = open(psmiclean, "w")
                    fsmiclean.write("ERROR")
                    fsmiclean.close()
                    return 1

                else:
                    self.log = self.log + "Prepared SMI :" + str(smilesclean) + "\n"

                    self.smiclean = smilesclean

                    # convert in inchikey
                    molformat = Chem.MolFromSmiles(smilesclean)
                    inchi = Chem.inchi.MolToInchi(molformat)
                    inchikey = Chem.inchi.InchiToInchiKey(inchi)
                    self.inchikey = inchikey
                    self.psmiclean = psmiclean
                    self.smiclean = smilesclean

                    fsmiclean = open(psmiclean, "w")
                    if inchikey != "" :
                        fsmiclean.write(smilesclean + "\t" + inchikey + "\t" + self.name)
                        fsmiclean.close()
                        return 0
                    else:
                        fsmiclean.write("ERROR")
                        self.err = 1
                        fsmiclean.close()
                        return 1



    def compute1D2DDesc(self, prDescbyChem):

        self.prDesc1D2D = prDescbyChem
        # check if descriptors already computed
        pdes = prDescbyChem + self.inchikey + ".txt"
        if self.inchikey == "":
            self.err = 1
            self.log = self.log + "No clean SMILES\n"
            return 1

        if path.exists(pdes) and path.getsize(pdes) > 10:
            filin = open(pdes, "r")
            llines = filin.readlines()
            filin.close()
            ldesc = llines[0].strip().split("\t")[1:]
            lval = llines[1].strip().split("\t")[1:]
            ddes = {}
            i = 0
            while i < len(ldesc):
                ddes[ldesc[i]] = lval[i]
                i += 1
            self.allDesc = ddes
            self.log = self.log + "Desc already computed -> " + pdes + "\n"
            return

        if not "smiclean" in self.__dict__:
            self.log = self.log + "No smiles prepared\n"
            self.err = 1
        else:
            self.mol = loader.ReadMolFromSmile(self.smiclean)

            try:
                self.consti = constitution.GetConstitutional(self.mol)
            except:
                self.consti = {}
            self.compo = {}
            try:
                self.compo["nheavy"] = self.mol.GetNumHeavyAtoms()
            except:
                self.compo = {}

            try:
                self.molprop = molproperty.GetMolecularProperty(self.mol)
            except:
                self.molprop = {}

                # 2D
            try:
                self.topo = topology.GetTopology(self.mol)
            except:
                self.topo = {}
            try:
                self.connect = connectivity.GetConnectivity(self.mol)
            except:
                self.connect = {}
            try:
                self.kap = kappa.GetKappa(self.mol)
            except:
                self.kap = {}
            try:
                self.burden = bcut.GetBurden(self.mol)
            except:
                self.burden = {}
            try:
                self.basakD = basak.Getbasak(self.mol)
            except:
                self.basakD = {}
            try:
                self.est = estate.GetEstate(self.mol)
            except:
                self.est = {}
            try:
                self.moreauBurto = moreaubroto.GetMoreauBrotoAuto(self.mol)
            except:
                self.moreauBurto = {}
            try:
                self.autcormoran = moran.GetMoranAuto(self.mol)
            except:
                self.autcormoran = {}
            try:
                self.gearycor = geary.GetGearyAuto(self.mol)
            except:
                self.gearycor = {}
            try:
                self.charges = charge.GetCharge(self.mol)
            except:
                self.charges = {}
            try:
                self.MOE = moe.GetMOE(self.mol)
            except:
                self.MOE = {}

            # combine all 1D2D
            if not "allDesc1D2D" in dir(self):
                self.allDesc1D2D = dict()
            self.allDesc1D2D.update(deepcopy(self.consti))
            self.allDesc1D2D.update(deepcopy(self.compo))
            self.allDesc1D2D.update(deepcopy(self.molprop))
            self.allDesc1D2D.update(deepcopy(self.topo))
            self.allDesc1D2D.update(deepcopy(self.connect))
            self.allDesc1D2D.update(deepcopy(self.kap))
            self.allDesc1D2D.update(deepcopy(self.burden))
            self.allDesc1D2D.update(deepcopy(self.basakD))
            self.allDesc1D2D.update(deepcopy(self.est))
            self.allDesc1D2D.update(deepcopy(self.moreauBurto))
            self.allDesc1D2D.update(deepcopy(self.autcormoran))
            self.allDesc1D2D.update(deepcopy(self.gearycor))
            self.allDesc1D2D.update(deepcopy(self.charges))
            self.allDesc1D2D.update(deepcopy(self.MOE))

    #to write
    def compute3DDesc(self, pr3DDesc):

        self.prDesc3D = pr3DDesc
        p3Ddesc = pr3DDesc + self.inchikey + ".txt"
        if self.inchikey == "":
            return 1
        # control already compute
        if path.exists(p3Ddesc)and path.getsize(p3Ddesc) > 10:
            filin = open(p3Ddesc, "r")
            llines = filin.readlines()
            filin.close()
            ldesc = llines[0].strip().split("\t")[1:]
            lval = llines[1].strip().split("\t")[1:]
            ddes = {}
            i = 0
            while i < len(ldesc):
                ddes[ldesc[i]] = lval[i]
                i += 1
            self.allDesc3D = ddes
            self.log = self.log + "Desc already computed -> " + p3Ddesc + "\n"
            return 0

        else:
            if not "psdf3D" in self.__dict__:
                self.log = self.log + "ERROR no 3D chemical file\n"
                self.err = 1
                return 1

            ddesc = descriptors3D.get3Ddesc(self.psdf3D)
            if ddesc == {}:
                self.log = self.log + "ERROR computation 3D descriptors\n"
                return 1
            else:
                self.allDesc3D = ddesc
                return 0





    def writeDesc(self, ldesc, filin, ltypeDesc):

        lw = []
        flag = 0
        for desc in ldesc:
            for typeDesc in ltypeDesc:
                if typeDesc == "1D2D":
                    if "allDesc1D2D" in self.__dict__:
                        if desc in self.allDesc1D2D:
                            try: lw.append(str(self.allDesc1D2D[desc]))
                            except: lw.append("NA")
                            flag = 1
                elif typeDesc == "3D":
                    if "allDesc3D" in self.__dict__:
                        if desc in self.allDesc3D:
                            try: lw.append(str(float(self.allDesc3D[desc])))
                            except: lw.append("NA")
                            flag = 1
            if flag == 0:
                lw.append("NA")
            else:
                flag = 0

        filin.write(self.name + "\t" + "\t".join(lw) + "\n")



    def writelog (self, prout):

        plog = prout + self.name + ".log"
        flog = open(plog, "w")
        flog.write(self.log)
        flog.close()





    def writeTablesDesc(self, prDescbyCAS, typeDesc):

        if typeDesc == "1D2D":
            if "allDesc1D2D" in self.__dict__:

                ptable = prDescbyCAS + self.inchikey + ".txt"
                ftable = open(ptable, "w")
                ftable.write("ID\t" + "\t".join(self.allDesc1D2D.keys()) + "\n")
                ftable.write(self.inchikey)
                for desc in self.allDesc1D2D.keys():
                    if self.allDesc1D2D[desc] == "inf" or self.allDesc1D2D[desc] == "nan":
                        self.allDesc1D2D[desc] = "NA"
                    ftable.write("\t" + str(self.allDesc1D2D[desc]))
                ftable.write("\n")
                ftable.close()
                return 0
            else:
                self.log = self.log + "No descriptors computed 1D2D for table\n"
                return 1

        if typeDesc == "3D":
            if "allDesc3D" in self.__dict__:

                ptable = prDescbyCAS + self.inchikey + ".txt"
                ftable = open(ptable, "w")
                ftable.write("ID\t" + "\t".join(self.allDesc3D.keys()) + "\n")
                ftable.write(self.inchikey)
                for desc in self.allDesc3D.keys():
                    if self.allDesc3D[desc] == "inf" or self.allDesc3D[desc] == "nan":
                        self.allDesc3D[desc] = "NA"
                    ftable.write("\t" + str(self.allDesc3D[desc]))
                ftable.write("\n")
                ftable.close()
                return 0
            else:
                self.log = self.log + "No descriptors computed 1D2D for table\n"
                return 1

        # maybe add all


    def generate3DFromSMILES(self, prSDF3D, software):
        """
        Compute descriptors 3D from SMILES code and generate the 3D in sdf format save in prSDF3D
        :return: dictionary of descriptors in all3D
        """
        self.software3D = software

        #log
        self.log = self.log + "Generate 3D using -> " + software + "\n"

        # filout
        psdf3D = prSDF3D + self.name + ".sdf"

        # control if chemical exists
        #print psdf3Dout
        if path.exists(psdf3D):
            self.log = self.log + "3D already generated\n"
            self.psdf3D = psdf3D
            self.err = 0

        if software == "ligprep":

            # Create temporary folder with smi inside
            prtemp = pathFolder.createFolder(prSDF3D + "temp3D/", clean=1)
            pfilesmile = prtemp + "tem.smi"
            filesmile = open(pfilesmile, "w")
            filesmile.write(self.smiclean)
            filesmile.close()

            # run ligprep
            psdf3Dtemp = runExternalSoft.runLigprep(psmilin=pfilesmile)
            # case error in ligprep
            if not path.exists(psdf3Dtemp) or path.getsize(psdf3Dtemp) == 0:
                self.log = self.log + "EROOR in ligPrep generation\n"
                self.err=1
            else:
                toolbox.selectMinimalEnergyLigPrep(psdfin=psdf3Dtemp,
                                                               psdfout=psdf3D)
                # take only best energy
                pathFolder.cleanFolder(prtemp)


        elif software == "RDKit":
            #generation using the method of Riniker and Landrum

            if not "smiclean" in self.__dict__:
                self.log = self.log + "ERROR in rdkit 3D generation not clean smiles\n"
                self.err = 1
            else:
                mol = Chem.MolFromSmiles(self.smiclean)
                molH = Chem.AddHs(mol)
                err = AllChem.EmbedMolecule(molH, AllChem.ETKDG())
                if err == 1:
                    self.log = self.log + "ERROR in rdkit 3D generation\n"
                    self.err = 1
                else:
                    wmol = Chem.MolToMolBlock(molH)
                    pmol = psdf3D[0:-3] + "mol"
                    fmol3D = open(psdf3D[0:-3] + "mol", "w")
                    fmol3D.write(wmol)
                    fmol3D.close()

                    runExternalSoft.babelConvertMoltoSDF(pmol, psdf3D)

                    if path.exists(psdf3D) and path.getsize(psdf3D) > 100:
                        self.psdf3D = psdf3D
                        remove(pmol)
                    else:
                        self.log = self.log + "ERROR during the write process of sdf\n"
                        self.err = 1

######## not Use



    def generate3DLigPrepFromSMILES(self, log):
        """
        Compute descriptors 3D from SMILES code and generate the 3D using ligprep
        :return: dictionary of descriptors in all3D
        """

        # define folder with selected 3D structure
        pr3DSDF = pathFolder.createFolder(self.prout + "SDF3D/")

        # clean temp folder - used to compute 3D descriptors
        prtemp = pathFolder.createFolder(self.prout + "temp3D/", clean = 1)
        psdf3Dout = pr3DSDF + self.compound[self.namek] + ".sdf"

        # control if chemical exists
        #print psdf3Dout
        if path.exists(psdf3Dout):
            return pr3DSDF

        # temp SMILES
        pfilesmile = prtemp + "tem.smi"
        filesmile = open(pfilesmile, "w")
        filesmile.write(self.compound["SMILES"])
        filesmile.close()

        pdesc = ""
        # run ligprep
        if not path.exists(psdf3Dout):
            psdf3D = runExternalSoft.runLigprep(psmilin=pfilesmile)
            #ffff
            # case error in ligprep
            if not path.exists(psdf3D) or path.getsize(psdf3D) == 0:
                log.write(self.compound[self.namek] + "\t" + self.compound["SMILES"] + "\t" + psdf3D)
            else:
                psdf3Dout = toolbox.selectMinimalEnergyLigPrep(psdfin=psdf3D,
                                                               psdfout=psdf3Dout)
                # take only best energy
                pathFolder.cleanFolder(prtemp)
        return pr3DSDF







#### not sure we need 3D


#def get_descriptor3Down(pr3DSDF, pdesc3D, geo=1, cpsa=1, rdf=1, morse=1, whim=1):

#    i = 0
#    lsdf = listdir(pr3DSDF)
#    ddesc = {}
#    for psdf in lsdf[:7400]:
#        if psdf[-3:] == "sdf":
#            name = psdf[:-4]
#            print i, name, "3D"
#            ddesc[name]=descriptors3D.get3Ddesc(pr3DSDF + psdf, geometry=geo, cpsa=cpsa, rdf=rdf, morse=morse, whim=whim)
#        i += 1
#    lheader = ddesc[name].keys()
#
    # write table desc
#    filout = open(pdesc3D, "w")
#    filout.write("ID\t" + "\t".join(lheader) + "\n")

#    for name in ddesc.keys():
#        filout.write(str(name))
#        for desc in lheader:
#            filout.write("\t" + str(ddesc[name][desc]))
#        filout.write("\n")
#    filout.close()

#    return pdesc3D
