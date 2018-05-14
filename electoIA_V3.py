#!/usr/bin/env python2
# Name: Ioannis Anastopoulos
# Date: 04/18/2018

import numpy as np
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
from matplotlib.lines import Line2D
import time
import scipy
from scipy import stats
from random import shuffle
from random import sample
from multiprocessing import Pool

clinical_file ="/projects/sysbio/users/TCGA/PANCANATLAS/Data/Clinical/out.tsv"
id_mapping_file ="/projects/sysbio/users/TCGA/PANCANATLAS/Data/ID_mapping.tsv"
sim_table="/soe/vfriedl/share/Pancanatlas_TM_similarity_matrix.tsv"
mutation_file="/projects/sysbio/users/TCGA/PANCANATLAS/Data/SNVs/Ver1.3/out.tsv"

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
				print 'no samples in ranked_samples'
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
				distance = abs(bin_ - bin2)*bins

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
	def testSeparation(pos_samples, neg_samples,pseudocounts,pos_distribution,neg_distribution, bins):
		pos_infered_data = []
		neg_infered_data = []
		number_pos_samples = len(pos_samples) + pseudocounts
		number_neg_samples = len(neg_samples)

		breaks = np.linspace(0,1,bins+1)

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
		print 'negKS_distribution',negKS_distribution
		# ks_distance, ks_pvalue, direction = self.testSeparationRaw(posKS_distribution, negKS_distribution)
		bins=50
		pseudocounts=2
		dummy_prior=0.1
		neg_hist,neg_bins = self.makeHistograms(negKS_distribution, bins)
		print 'raw_histo',neg_hist
		# print('pos bins from np.histogram',pos_hist,min(pos_bins), max(pos_bins))
		neg_Smoothed_histo=self.smoothHistogram(neg_hist, bins) #smoothed counts

		return neg_Smoothed_histo, neg_bins
	
	def returnSmoothPosHistogram(self, neg_Smoothed_histo):
		posKS_distribution=self.KS_distribution(self.pos_samples)
		# ks_distance, ks_pvalue, direction = self.testSeparationRaw(posKS_distribution, negKS_distribution)
		bins=50
		pseudocounts=2
		dummy_prior=0.1
		pos_hist,pos_bins = self.makeHistograms(posKS_distribution, bins)
		# print('pos bins from np.histogram',pos_hist,min(pos_bins), max(pos_bins))
		pos_pseudo_histo=self.addPseudocounts(pos_hist,neg_Smoothed_histo,pseudocounts ) #adding pseudocounts to positive distribution
		pos_pseudoSmoothed_histo=self.smoothHistogram(pos_pseudo_histo, bins) #smoothing pseudocounted positive distribution

		return pos_pseudoSmoothed_histo, pos_bins

	



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

def testSeparation(neg_samples,new_negsamples,smoothNegDistribution,new_smoothNegDistribution, bins):
	neg_infered_data = []
	new_neg_infered_data = []
	number_neg_samples = len(neg_samples)
	number_new_neg_samples = len(new_negsamples)

	breaks = np.linspace(0,1,bins+1)

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

def negDownsampling(pos_samples,neg_samples, sim_table, feature):
	shuffle(neg_samples)
	# print 'starting downsampling', (time.time()-startTime)/60
	n_samples=100
	neg_splits=np.array_split(list(neg_samples), len(neg_samples)/n_samples)

	ks_2samp_old=None 
	ks_threshold=0.004 #0.01 0.001 0.0001
	
	first_split=list(neg_splits[0])
	split1=first_split

	SamplesRanking= simMatrix(sim_table,pos_samples, split1) #ranking of pos and negative sampels for the feature we are testing
	electo1=Electo(SamplesRanking,list(pos_samples), split1)
	smoothNegDistribution, neg_bins=electo1.returnSmoothNegHistogram()
	

	ks_distances=[]
	ks_p=[]
	smoothNeg=smoothNegDistribution
	total_iterations=0
	with open('./test_logs/'+feature+'_Downsampling_test_with_%d_negSamples.txt' %n_samples,'w') as f:
		f.write('\t'+'Old split'+'\t'+'new split lenght'+'\t'+'2samp ks'+'\n')

		for i in range(1,len(neg_splits)):	
			new_split=split1+list(neg_splits[i]) #next chunk
			
			new_SamplesRanking= simMatrix(sim_table,pos_samples, list(new_split))
			electo2=Electo(new_SamplesRanking,list(pos_samples), list(new_split))
			new_smoothNegDistribution, new_neg_bins=electo2.returnSmoothNegHistogram()

			ks=testSeparation(split1,new_split,smoothNeg,new_smoothNegDistribution,50)
			print 'Iteration %d split lenght' %i,len(split1),'new split lenght',len(new_split),'2samp ks',ks
			
			ks_distances.append(ks[0])
			ks_p.append(ks[1])
			total_iterations+=1
			f.write('Iteration {} split length \t {} \t {} \t {} \n'.format(i,len(split1),len(new_split),ks))
			f.write('time to downsample \t {}mins \n'.format((time.time()-startTime)/60))
			
			ks_delta=+float('inf')
			if ks_2samp_old is not None:
				ks_delta=abs(ks_2samp_old-ks[0])
			
			if ks[0] <0.05 and ks_delta<=ks_threshold:
				smoothNeg=new_smoothNegDistribution
				new_neg_bins=new_neg_bins #this may be redundant 
				new_ranking=new_SamplesRanking
				break
			else:
				ks_2samp_old=ks[0]
				split1=new_split
				smoothNeg=new_smoothNegDistribution
				continue
			
		f.write('ks delta used {} \n'.format(ks_threshold))
		f.write('total iterations {} \n'.format(total_iterations))
	# return ks_distances, ks_p
	return smoothNeg, new_split,new_neg_bins,new_ranking


def preProcessingDownSampling():
	starttime=time.time()
	print 'im working',time.time()-starttime

	#feature = 'RNF14'
	for feature in ['NRAS']:
		fig,ax = plt.subplots(1)
		feature_dict, all_patients = file_parser(clinical_file, mutation_file, feature)
		pos_samples, neg_samples=getPosNegsamples(feature, feature_dict, all_patients) #pos and neg samples for each feature
		sampleRanks=simMatrix(sim_table, list(pos_samples), list(neg_samples))
		#IMPORTANT: get intersection of samples between clinical, mutation, and sim_table to avoid KeyWord errors
		pos_samples=set(pos_samples).intersection(sampleRanks.keys())
		neg_samples=set(neg_samples).intersection(sampleRanks.keys())
		print feature, 'pos samples',len(pos_samples), 'neg samples', len(neg_samples)
		# print 'pos_neg and all ranking done', (time.time()-starttime)/60, 'min'

		electo=Electo(sampleRanks,list(pos_samples), list(neg_samples)) #3 minutes faster with just 100 samples 
		all_neg_histo, all_neg_bins= electo.returnSmoothNegHistogram()
		smoothed_pos, pos_bins=electo.returnSmoothPosHistogram(all_neg_histo)
		ks_d, ks_p=electo.testSeparation(list(pos_samples), list(neg_samples),2,smoothed_pos,all_neg_histo,50)

		shuffle(list(neg_samples))
		# print('neg_samples', len(neg_samples))
		chunk_sizes=np.arange(10,1010,10)
		plt.plot([0, 1005], [ks_d, ks_d], color='black',linewidth=1.1, zorder=3)
		plt.text(max(chunk_sizes)+5, ks_d+0.01, 'Pos-Neg distance', fontsize=6, horizontalalignment='right')
		for n_samples in chunk_sizes:
			neg_splits=np.array_split(list(neg_samples), len(neg_samples)/n_samples)

			ks_pos=[] #ks distances between the chunk and whole positive
			ks_neg=[] #ks distances betweent eh chunk and whole negative
			
			for chunk in neg_splits:
				print chunk
				electo_chunk=Electo(sampleRanks,list(pos_samples), list(chunk)) #3 minutes faster with just 100 samples 
				chunk_neg_histo, chunk_neg_bins= electo_chunk.returnSmoothNegHistogram()
				print 'NRAS neg distribution', chunk_neg_histo
				ks_d1, ks_p1=electo.testSeparation(list(pos_samples), list(chunk),2,smoothed_pos,chunk_neg_histo,50)
				ks_pos.append(ks_d1)
				ks_d2, ks_p2=testSeparation(list(neg_samples),list(chunk),all_neg_histo,chunk_neg_histo,50)
				ks_neg.append(ks_d2)
			ks_pos_mean=np.mean(ks_pos)
			ks_pos_sd=np.std(ks_pos)
			ks_neg_mean=np.mean(ks_neg)
			ks_neg_sd=np.std(ks_neg)

			# ax.scatter(n_samples,ks_pos_sd,s=3.5,c='blue',edgecolors='blue', marker='v',linewidth=0,label='Pos p-val')
			# ax.scatter(n_samples,ks_neg_sd,s=3.5,c='black',edgecolors='black', marker='v',linewidth=0,label='Neg p-val')

			ax.errorbar(n_samples, ks_pos_mean, yerr=ks_pos_sd, ms=2.5, linewidth=0, mfc='blue',mec='blue', marker='o',
									ecolor='red', elinewidth=1, capsize=1.5, label='KS distance to positive distribution')
			ax.errorbar(n_samples, ks_neg_mean, yerr=ks_neg_sd, ms=2.5, linewidth=0, mfc='black',mec='black', marker='o',
								ecolor='green', elinewidth=1, capsize=1.5, label='KS distance to negative distribution')
			# print 'plotted chunk'+str(n_samples), (time.time()-starttime)/60, 'min'

		handles, labels = plt.gca().get_legend_handles_labels()
		by_label = OrderedDict(zip(labels, handles))
		ax.legend(by_label.values(), by_label.keys(),loc='upper right')
		plt.xticks(list(np.arange(10,1010,50)),list(map(lambda x:str(x),list(np.arange(10,1010,50)))), fontsize=6)
		
		ax.set_xlim(0,max(chunk_sizes)+5)
		ax.set_ylabel(r'KS distance')
		ax.set_xlabel(r'Negative Chunk Size')
		plt.title(feature+' KS distances by chunk size')
		plt.savefig('./plots/'+feature+'PreprocessDownsampling.png', dpi=600)
		print 'DONE with'+' feature', (time.time()-starttime)/60, 'min'



def chunk(feature):
	global startTime
	startTime = time.time()
	print 'starting downsampling and 2_samp ks test', time.time()-startTime
	# feature='RNF14'
	feature_dict, all_patients = file_parser(clinical_file, mutation_file, feature)
	pos_samples, neg_samples=getPosNegsamples(feature, feature_dict, all_patients) #pos and neg samples for each feature
	new_neg_histo,new_neg_samples, new_neg_bins,new_ranking=negDownsampling(list(pos_samples),list(neg_samples), sim_table, feature)
	electo=Electo(new_ranking,list(pos_samples), list(new_neg_samples))
	smoothPosHisto, pos_bins=electo.returnSmoothPosHistogram(new_neg_histo)

	ks_d1, ks_p1=electo.testSeparation(list(pos_samples), list(new_neg_samples),2,smoothPosHisto,new_neg_histo,50)
	print ks_d1, ks_p1
	print 'DONE' ,(time.time()-startTime)/60,'min'
	return ks_d1, ks_p1, (time.time()-startTime)/60, new_neg_samples,new_neg_histo
def normal(feature, size=None):

	global startTime_allneg
	startTime_allneg=time.time()
	print 'starting normal 2_samp ks test', time.time()-startTime_allneg
	# feature='RNF14'
	feature_dict, all_patients = file_parser(clinical_file, mutation_file, feature)
	
	pos_samples, neg_samples=getPosNegsamples(feature, feature_dict, all_patients) #pos and neg samples for each feature
	if size is not None:
		neg_samples=sample(list(neg_samples),int(size))
	print 'testing for neg size', len(neg_samples)
	sampleRanks=simMatrix(sim_table, list(pos_samples), list(neg_samples))
	

	electo=Electo(sampleRanks,list(pos_samples), list(neg_samples))
	smoothNegHisto,neg_bins=electo.returnSmoothNegHistogram()
	smoothPosHisto, pos_bins=electo.returnSmoothPosHistogram(smoothNegHisto)
	

	ks_d1, ks_p1=electo.testSeparation(list(pos_samples), list(neg_samples),2,smoothPosHisto,smoothNegHisto,50)
	print ks_d1, ks_p1
	print 'DONE' ,(time.time()-startTime_allneg)/60,'min'
	return ks_d1, ks_p1, (time.time()-startTime_allneg)/60,neg_samples,smoothNegHisto

def chunkVSnormal():
	fig,ax = plt.subplots(1)
	custom_lines = [Line2D([0], [0], color='yellow', lw=2),
					Line2D([0], [0], color='red', lw=2),
					Line2D([0], [0], color='green',lw=2),
					Line2D([0], [0], color='blue',lw=2)]

	def bar_plot(hist, bins, color='black'):
		for values in np.arange(0,len(hist)):
			left=bins[values]
			bottom=0
			height=hist[values]
			width=bins[values+1]-bins[values]
			rectangle=mplpatches.Rectangle((left,bottom),width,height,
											facecolor=color,
											edgecolor='black',
											linewidth=0.9,
											alpha=0.2) 
			ax.add_patch(rectangle)

		
	features= ['RNF14','BRAF','IDH1','VHL','KRAS','NRAS']
	with open('./test_logs/Chunk_vs_Normal.txt','w') as f:
		f.write('Feature \t Neg_Chunk_Size \t Chunk_Downsampling_Time(min) \t Chunk_ks_distance \t Chunk_p_val \t All_neg_size \t Normal_time(min) \t Normal_ks_distance \t Normal_p_val \t 2samp_d:Chunk_vs_Normal \t 2samp_p_val:Chunk_vs_Normal \n')
		for feature in features:

			ks_d_chunk, ks_p_chunk, chunk_time,Neg_Chunk_Samples,new_neg_histo,new_neg_bins, new_pos_histo,new_pos_bins=chunk(feature)

			ks_d, ks_p, time, neg_samples,smoothNegHisto,neg_bins,pos_histo,pos_bins=normal(feature)
			two_samp_ks_d, two_samp_ks_p=testSeparation(neg_samples, Neg_Chunk_Samples,smoothNegHisto,new_neg_histo, 50)
			f.write('{} \t {} \t {:.3f} \t {:.3f} \t {:.3f} \t {} \t {:.3f} \t {:.3f} \t {:.3f} \t {:.3f} \t {:.3f} \n'.format(feature,len(Neg_Chunk_Samples),chunk_time,ks_d_chunk, ks_p_chunk,  len(neg_samples),time,ks_d, ks_p ,two_samp_ks_d,two_samp_ks_p))

			bar_plot(smoothNegHisto,neg_bins, color='red')
			bar_plot(new_neg_histo,new_neg_bins, color='yellow')
			bar_plot(pos_histo,pos_bins, color='green')
			bar_plot(new_pos_histo,new_pos_bins, color='blue')

			ax.legend(custom_lines, ['Chunked Neg Distribution', 'Whole Neg Distribution', 'Pos Distribution from all_neg','Pos Distribution from chunk_neg'])	
			ax.set_ylim(0,np.max([np.max(smoothNegHisto),np.max(new_neg_histo),np.max(pos_histo),np.max(new_pos_histo)])+0.05)
			plt.title(feature+'_WholeNegDist_vs_ChunkedNegDist')
			plt.savefig('./plots/%s_WholeNegDist_vs_ChunkedNegDist.png' %feature, dpi=600)

def preprocesingVSnromal():
	
	features= [('KRAS',320),('IDH1',200),('RNF14',320),('BRAF',260)]
	with open('./test_logs/Preprocessing_vs_Normal.txt','w') as f:
		f.write('Feature \t Preprocessing_Size \t Preprocessing_Time(min) \t Preprocessing_ks_distance \t Preprocessing_p_val \t All_neg_size \t Normal_time(min) \t Normal_ks_distance \t Normal_p_val \t 2samp_d:Chunk_vs_Normal \t 2samp_p_val:Chunk_vs_Normal \n')
		for feature in features:

			ft = feature[0]
			n_neg_samples = feature[1]

			ks_d_chunk, ks_p_chunk, chunk_time,Neg_Chunk_Samples,new_neg_histo=normal(ft,n_neg_samples)

			ks_d, ks_p, time, neg_samples,smoothNegHisto=normal(ft)
			two_samp_ks_d, two_samp_ks_p=testSeparation(neg_samples, Neg_Chunk_Samples,smoothNegHisto,new_neg_histo, 50)
			f.write('{} \t {} \t {:.3f} \t {:.3f} \t {:.3f} \t {} \t {:.3f} \t {:.3f} \t {:.3f} \t {:.3f} \t {:.3f} \n'.format(feature,len(Neg_Chunk_Samples),chunk_time,ks_d_chunk, ks_p_chunk,  len(neg_samples),time,ks_d, ks_p ,two_samp_ks_d,two_samp_ks_p))
			
def main():
	preProcessingDownSampling()

if __name__ == "__main__":
	main()


# Prediction of features by similarity - Pro.Fe.Si.








