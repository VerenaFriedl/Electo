#!/usr/bin/env python2
# Name: Ioannis Anastopoulos
# Date: 04/18/2018

import numpy as np
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import time
import scipy
from scipy import stats
from random import shuffle
from multiprocessing import Pool

class Electo():

	def __init__(self, ranked_samples, pos_samples, neg_samples):
		self.ranked_samples=ranked_samples #ranked samples based on their similarity to each other sample in the cohort by Spearman Correlation
		self.pos_samples=pos_samples #positive samples for the feature we are examining
		self.neg_samples=neg_samples #negative samples for the feature we are examining
	
	def KS_distribution(self,posORnegSamples):
		ks_distances=[]
		for sample in posORnegSamples:
			rank_list=[]
			#rank of samples based on positive samples
			if sample in self.ranked_samples: #this takes care of the case where the mRNA ID for the sample is not in the 
											  #similarity file, but it is in the id_mapping file, which is where the feature
											  #dictionary and mRNA SNV IDs are collected
				ranked_samples=list(self.ranked_samples[sample])
				# print sample, ranked_samples
				for s in self.pos_samples:
					if s in ranked_samples and s!=sample:
						rank_list.append(ranked_samples.index(s))
				
				ks = scipy.stats.kstest(rank_list, 'uniform',args=(0,len(ranked_samples)-1),alternative="greater")
				ks_distances.append(ks[0])
			else:
				continue
		
		return ks_distances

	def testSeparationRaw(self,posKS_distribution,negKS_distribution):
		""" 2 sample KS test between the positive and negative raw distributions"""

		ks = scipy.stats.ks_2samp(posKS_distribution,negKS_distribution)
		ks_distance = ks[0]
		ks_pvalue = ks[1]
		direction = 1
		return ks_distance, ks_pvalue, direction
		
	def makeHistograms(self,KS_distribution,bins):
		"""makes histograms from distributions"""

		hist,hist_bins=np.histogram(KS_distribution,np.linspace(0,1,bins+1))

		return hist,hist_bins
	def smoothHistogram(self,counts, bins):
		"""smooth a histogram
		input: the histogram counts, smoothing parameter alpha (defaults to 1)
		return: smoothed histogram counts (normalized to sum up to 1)"""

		counts_smoothed = []
		alpha = 1
		breaks = np.linspace(0,1,bins+1)

		for i in range(1,bins+1):
			bin_ = breaks[i]
			bin_sum = 0
			for j in range(1, bins+1):
				count = counts[j-1]
				bin2 = breaks[j]
				distance = abs(bin_ - bin2)*20

				bin_sum += (2**(-alpha*distance))*count

			counts_smoothed.append(bin_sum)
		counts_smoothed_= list(map(lambda x:x/sum(counts_smoothed), counts_smoothed))
		return counts_smoothed_ 
	
	def addPseudocounts(self,pos_histogram, neg_smoothed, pseudocounts):
		""" Add pseudocounts to histogram bins"""
		counts_pseudo =[]
		for i in range(len(neg_smoothed)):
			base_frequency = neg_smoothed[i]
			new_count = float(pos_histogram[i]) + (float(base_frequency) * int(pseudocounts))
			counts_pseudo.append(new_count)

		return counts_pseudo

	@staticmethod
	def testSeparation(pos_samples, neg_samples,pseudocounts,pos_distribution,neg_distribution):
		pos_infered_data = []
		neg_infered_data = []
		number_pos_samples = len(pos_samples) + pseudocounts
		number_neg_samples = len(neg_samples)

		breaks = np.linspace(0,1,21)

		for i in range(len(pos_distribution)):
			pos_count=int(round(pos_distribution[i]*number_pos_samples))
			neg_count=int(round(neg_distribution[i]*number_neg_samples))

			pos_random = np.random.uniform(low=breaks[i], high=breaks[i+1], size=pos_count)
			neg_random = np.random.uniform(low=breaks[i], high=breaks[i+1], size=neg_count)

			pos_infered_data+=(pos_random.tolist())			
			neg_infered_data+=(neg_random.tolist())
		ks = scipy.stats.ks_2samp(pos_infered_data,neg_infered_data)
		ks_distace=ks[0]
		ks_p=ks[1]
		return ks_distace, ks_p
	@staticmethod
	def LogOddsRatio(neg_distribution,pos_distribution):
		log_ratios = []

		for i in range(1,len(neg_distribution)):
			p_neg = neg_distribution[i]
			p_pos = pos_distribution[i]
			log_ratios.append(np.log10(p_pos/p_neg))

		return log_ratios
	@staticmethod
	def Probabilities(prior,neg_distribution,pos_distribution):
		probability_vector = []
		for i in range(len(neg_distribution)):
			p_neg = neg_distribution[i]
			p_pos = pos_distribution[i]
			prob = (p_pos*prior) / (p_pos*prior + (p_neg*(1-prior)))
			probability_vector.append(prob)

		return probability_vector

	def returnSmoothNegHistogram(self):
		negKS_distribution=self.KS_distribution(self.neg_samples)
		# ks_distance, ks_pvalue, direction = self.testSeparationRaw(posKS_distribution, negKS_distribution)
		bins=20
		pseudocounts=2
		dummy_prior=0.1
		neg_hist,neg_bins = self.makeHistograms(negKS_distribution, bins)
		# print('pos bins from np.histogram',pos_hist,min(pos_bins), max(pos_bins))
		neg_Smoothed_histo=self.smoothHistogram(neg_hist, bins) #smoothed counts

		return neg_Smoothed_histo,neg_bins
	
	def returnSmoothPosHistogram(self, neg_Smoothed_histo):
		posKS_distribution=self.KS_distribution(self.pos_samples)
		# ks_distance, ks_pvalue, direction = self.testSeparationRaw(posKS_distribution, negKS_distribution)
		bins=20
		pseudocounts=2
		dummy_prior=0.1
		pos_hist,pos_bins = self.makeHistograms(posKS_distribution, bins)
		# print('pos bins from np.histogram',pos_hist,min(pos_bins), max(pos_bins))
		pos_pseudo_histo=self.addPseudocounts(pos_hist,neg_Smoothed_histo,pseudocounts ) #adding pseudocounts to positive distribution
		pos_pseudoSmoothed_histo=self.smoothHistogram(pos_pseudo_histo, bins) #smoothing pseudocounted positive distribution

		return pos_pseudoSmoothed_histo,pos_bins

	



def file_parser(clinical_file, feature_file, feature):
	'''Function takes in files and disease string we're interested in and returns dict of lists of ranked samples according to similarity with each other
	samples in the cohort, and a dict of features and the samples positive for that feature'''
	###############-------------getting patient IDs-----------############
	all_patients=set() #all patient IDs from clinical file
	with open(clinical_file) as f:
		lines=f.readlines()
		
		for line in lines[1:]:
			s=line.split()
			all_patients.add(s[1]) 
	###############-------------getting patient IDs------------############

	###############-------------getting mRNA and SNV IDs------------############
	# SNV_mRNA_dict={} #SNV ID to mRNA id? or are these ID patient IDs??
	# mRNA_ID=set() #mRNA IDs of all patients
	# SNV_ID=set() #SNV IDs of all aptients
	# with open(id_mapping_file) as f:
	# 	lines=f.readlines()
	# 	for line in lines[1:]:
	# 		line=line.split()
	# 		if line[0] in all_patients and line[3]!='NA' and line[6]!='NA': #Excluding NA from mRNA and SNV IDs
	# 			mRNA_ID.add(line[6])
	# 			SNV_ID.add(line[3])
	# 			SNV_mRNA_dict[line[3]]=line[6] #SNV id: mRNA id dict to be used in feature creation
	###############-------------getting mRNA and SNV IDs ------------############


	###############-------------getting positive samples for each feature------------------------############
	#should we make this an option?
	variant_type = ('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site')
	feature_lines=[] #list of lists
	feature_pos_samples= defaultdict(set)#{gene: set of mRNA_IDs positive for feature }
	with open(feature_file) as f:
		lines=f.readlines()
		for line in lines:
			line=line.split()
			if line[0]==feature and line[8] in variant_type  and str(line[15])[:12] in all_patients :
				feature_pos_samples[line[0]].add(str(line[15])[:12])#ONLY POSITIVE


	###############-------------getting positive samples for each feature------------------------############

	return feature_pos_samples, all_patients



	###############-------------creating comaprison matrix of mRNA IDs from sim_table------------############

def getPosNegsamples(feature, feature_dict, all_patients):
	'''function takes in feature and a feature dict, and returns a set of positive and negative samples for each feature'''
	pos_samples=feature_dict[feature]
	neg_samples=set(all_patients).difference(pos_samples)
	return pos_samples, neg_samples

def testSeparation(neg_samples,new_negsamples,smoothNegDistribution,new_smoothNegDistribution):
	neg_infered_data = []
	new_neg_infered_data = []
	number_neg_samples = len(neg_samples)
	number_new_neg_samples = len(new_negsamples)

	breaks = np.linspace(0,1,21)

	for i in range(len(smoothNegDistribution)):
		neg_count=int(round(smoothNegDistribution[i]*number_neg_samples))
		new_neg_count=int(round(new_smoothNegDistribution[i]*number_new_neg_samples))

		neg_random = np.random.uniform(low=breaks[i], high=breaks[i+1], size=neg_count)
		new_neg_random = np.random.uniform(low=breaks[i], high=breaks[i+1], size=new_neg_count)

		neg_infered_data+=(neg_random.tolist())			
		new_neg_infered_data+=(new_neg_random.tolist())

	ks = scipy.stats.ks_2samp(neg_infered_data,new_neg_infered_data)
	ks_distace=ks[0]
	ks_p=ks[1]
	return ks_distace, ks_p
	
def simMatrix(sim_table, pos_samples, neg_samples):
	all_samples=list(pos_samples)+list(neg_samples)
	SamplesRanking={}
	with open(sim_table) as f:
		lines=f.readlines()
		patient_IDs=list(lines[0].split()[1:]) #IDs in columns - should be same as IDs in rows since it is a similarity matrix
		for line in lines:
			patient_ID=line.split()[0][:12] #ID in each row
			if patient_ID in all_samples:
				SamplesRanking[patient_ID]=zip(patient_IDs,line.split()[1:])
				SamplesRanking[patient_ID]=sorted(SamplesRanking[patient_ID],key=lambda x:x[1], reverse=True)
				SamplesRanking[patient_ID] = [str(x[0][:12]) for x in SamplesRanking[patient_ID]]
	return SamplesRanking

def negDownsampling(pos_samples,neg_samples, sim_table):
	shuffle(neg_samples)
	print('starting downsampling', (time.time()-startTime)/60)
	n_samples=100
	neg_splits=np.array_split(list(neg_samples), len(neg_samples)/n_samples)

	ks_2samp_old=None 
	ks_threshold=0.004
	
	first_split=list(neg_splits[0])
	split1=first_split

	SamplesRanking= simMatrix(sim_table,pos_samples, split1) #ranking of pos and negative sampels for the feature we are testing
	electo=Electo(SamplesRanking,list(pos_samples), split1)
	smoothNegDistribution, neg_bins=electo.returnSmoothNegHistogram()
	

	ks_distances=[]
	ks_p=[]
	smoothNeg=smoothNegDistribution
	total_iterations=0
	with open('Downsampling_test_10_iterations_%d_samples' %n_samples,'w') as f:
		f.write('\t'+'Old split'+'\t'+'new split lenght'+'\t'+'2samp ks'+'\n')
		for i in range(1,len(neg_splits)):
			if i==21:
				break
			

			new_split=split1+list(neg_splits[i]) #next chunk
			
			new_SamplesRanking= simMatrix(sim_table,pos_samples, list(new_split))
			electo=Electo(new_SamplesRanking,list(pos_samples), list(new_split))
			new_smoothNegDistribution, new_neg_bins=electo.returnSmoothNegHistogram()

			ks=testSeparation(split1,new_split,smoothNeg,new_smoothNegDistribution)
			print('Iteration %d split lenght' %i,len(split1),'new split lenght',len(new_split),'2samp ks',ks)
			ks_distances.append(ks[0])
			ks_p.append(ks[1])
			
			f.write('Iteration {} split length \t {} \t {} \t {} \n'.format(i,len(split1),len(new_split),ks))
			f.write('time to downsample \t {}mins \n'.format((time.time()-startTime)/60))
			
			ks_delta=+float('inf')
			if ks_2samp_old is not None:
				ks_delta=abs(ks_2samp_old-ks[0])
			
			# if ks[0] <0.05 and ks_delta<=ks_threshold:
			# 	break
			# else:
			total_iterations+=1
			ks_2samp_old=ks[0]
			split1=new_split
			smoothNeg=new_smoothNegDistribution
				# continue
			
		# f.write('ks delta used {} \n'.format(ks_threshold))
		f.write('total iterations {} \n'.format(total_iterations))
	return ks_distances, ks_p
			


def main():
	global startTime
	startTime = time.time()
	print('im working...', time.time()-startTime)
	#Files
	clinical_file ="/projects/sysbio/users/TCGA/PANCANATLAS/Data/Clinical/out.tsv"
	id_mapping_file ="/projects/sysbio/users/TCGA/PANCANATLAS/Data/ID_mapping.tsv"
	sim_table="/soe/vfriedl/share/Pancanatlas_TM_similarity_matrix.tsv"
	mutation_file="/projects/sysbio/users/TCGA/PANCANATLAS/Data/SNVs/current/out.tsv"
	
	feature_dict, all_patients = file_parser(clinical_file, mutation_file, 'BRAF')
	pos_samples, neg_samples=getPosNegsamples('BRAF', feature_dict, all_patients) #pos and neg samples for each feature
	
	print('pos_nef samples done', (time.time()-startTime)/60, 'min')
	
	# sampleRanks=simMatrix(sim_table, list(pos_samples), list(neg_samples))
	# print sampleRanks[sampleRanks.keys()[0]][::-1]
	# print sampleRanks.keys()[0]
	# # print list(all_patients)[:4]
	# # for feature in feature_dict.keys():
	print len(pos_samples)
	print len(neg_samples)
	# after downsampling

	ks_distances, ks_p= negDownsampling(list(pos_samples),list(neg_samples), sim_table)

	iterations=list(range(1,21))
	plt.scatter(iterations,ks_distances, color='blue', label='ks_distances')
	plt.scatter(iterations,ks_p, color='black', label='p_values')
	
	plt.legend(bbox_to_anchor=(1, 0.15), loc=1, borderaxespad=0.)
	plt.savefig('./plots/test.png', dpi=600)
	# print('downsampling done', (time.time()-startTime)/60, 'min')
	# for k in (negDownsampling(list(pos_samples),list(neg_samples), sim_table)):
	# 	smoothNegDistribution=list(k[0])
	# 	neg_bins=list(k[1])
	# 	new_smoothNegDistribution=list(k[2])
	# 	new_neg_bins=list(k[3])
	# 	ks=k[4] 

	# 	print(smoothNegDistribution, neg_bins)
	


	# if len(pos_samples)!=1: IMPORTANT: SOME FEATURES HAVE NO POSITIVES
	# electo=Electo(sampleRanks,list(pos_samples), list(neg_samples)) #3 minutes faster with just 100 samples 
	# neg_Smoothed_histo, neg_bins= electo.returnSmoothNegHistogram()
	# pos_pseudoSmoothed_histo,pos_bins=electo.returnSmoothPosHistogram(neg_Smoothed_histo)
	# ks_distance, ks_p=electo.testSeparation(pos_samples, neg_samples,2,pos_pseudoSmoothed_histo,neg_Smoothed_histo)
	# print ks_distance
	# plotHistograms(pos_pseudoSmoothed_histo, neg_Smoothed_histo)
	# # else:
	# # 	continue
	print 'DONE', (time.time()-startTime)/60, 'min'
if __name__ == "__main__":
	main()
