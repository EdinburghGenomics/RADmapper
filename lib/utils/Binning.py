'''
Created on 25 Feb 2010

@author: tcezard
'''
import math
from utils.GenomeLoader import GenomeLoader
from utils import utils_logging
from utils import pass_funct
import logging
import operator


def frange(limit1, limit2 = None, increment = 1.):
    """
    Range function that accepts floats (and integers).
    
    Usage:
    frange(-2, 2, 0.1)
    frange(10)
    frange(10, increment = 0.5)
    
    The returned value is an iterator.  Use list(frange) for a list.
    """
    
    if limit2 is None:
        limit2, limit1 = limit1, 0.
    else:
        limit1 = float(limit1)
    count = int(math.ceil( (limit2 - limit1)/increment))
    return (limit1 + n*increment for n in xrange(count))

class Distribution_holder(object):
    """Hold the elements of a discrete distribution.
    efficient for distribution containing large number of element of the same value."""
    def __init__(self, update=False):
        '''
        Constructor
        '''
        self.distribution={}
        self.sum=0
        #self.sum2=0
        self.nb_elmt=0
        if update:
            self.S=0
            self.M=0
            self.add_value=self._add_value_update
            self.get_nb_element=self._get_nb_element_update
            self.get_sum=self._get_sum_update
            self.get_std_dev=self._get_std_dev_update
        else:
            self.add_value=self._add_value_no_update
            self.get_nb_element=self._get_nb_element_no_update
            self.get_sum=self._get_sum_no_update
            self.get_std_dev=self.get_real_std_dev
            
            
    def _add_value_no_update(self, value, weight=1):
        """Add one (or weight) value(s) to the distribution"""
        if self.distribution.has_key(value):
            self.distribution[value]+=weight
        else:
            self.distribution[value]=weight  
        
    def _add_value_update(self,value, weight=1):
        """Add one (or weight) value(s) to the distribution"""
        if self.distribution.has_key(value):
            self.distribution[value]+=weight
        else:
            self.distribution[value]=weight
        
        self.nb_elmt+=weight
        self.sum+=value*weight
        #self.sum2+=(value**2)*weight
        self._update_std_dev_add(value, weight)
        
        
    def remove_value(self,value, weight=1):
        if self.distribution.has_key(value) and self.distribution.get(value)>=weight:
            tmp_weight=self.distribution.get(value)
            tmp_weight-=weight
            if tmp_weight>0:
                self.distribution[value]=tmp_weight
            else:
                self.distribution.pop(value)
        else:
            if not self.distribution.has_key(value):
                raise Exception("Can't remove %s time %s from distribution key does not exist"%(weight,value))
            elif self.distribution.get(value)>=weight:
                raise Exception("Can't remove %s time %s from distribution only %s value in the distribution"%(weight,value,self.distribution.get(value)))
            else:
                raise Exception("It's really screwed up")
        self.nb_elmt-=weight
        self.sum-=value*weight
        #self.sum2-=(value**2)*weight
        self._update_std_dev_remove(value,weight)
        
            
    def get_binned_value(self, nb_bin):
        import numpy
        all_values,weights=self.get_sorted_value_and_weight()
        nb_per_bin, bin_edges=numpy.histogram(a=all_values, weights=weights, bins=nb_bin)
        bin_start = bin_edges[0]
        bins=[]
        for bin_end in bin_edges[1:]:
            bins.append(bin_start+(bin_end-bin_start)/2)
            bin_start=bin_end
        return bins,nb_per_bin
        
    def get_values(self):
        return self.distribution.keys()
    
    def get_sorted_value_and_weight(self, reverse=False):
        values=[]
        weights=[]
        all_tuples = sorted(self.distribution.iteritems(), key=operator.itemgetter(0),reverse=reverse)
        for value,weight in all_tuples:
            values.append(value)
            weights.append(weight)
        return values,weights
    
    def get_sorted_weight_and_value(self,reverse=False):
        values=[]
        weights=[]
        all_tuples =  sorted(self.distribution.iteritems(), key=operator.itemgetter(1), reverse=reverse)
        for value,weight in all_tuples:
            values.append(value)
            weights.append(weight)
        return values,weights
        
    def _get_nb_element_update(self):
        return self.nb_elmt
    
    def _get_nb_element_no_update(self):
        nb_element=0
        for value in self.distribution.keys():
            nb_element+=self.distribution.get(value)
            
        return nb_element
    
    def _get_sum_update(self):
        return self.sum
    
    def _get_sum_no_update(self):
        sum=0
        for value in self.distribution.keys():
            sum+=self.distribution.get(value)*value
        return sum
    
    def get_nb_of_value(self, value):
        return self.distribution.get(value)
    
    def get_distribution(self):
        return self.distribution
    
    
    
    def get_mean(self):
        """Calculate the mean assuming the distribution is normal."""
        nb_elmt=self.get_nb_element()
        sum=self.get_sum()
        if nb_elmt<1:
            return 0
        return sum/float(nb_elmt)
    
    def get_median(self):
        """Calculate the median."""
        results = self.get_percentiles([50])
        return results[0]
    
    def _update_std_dev_add(self, value ,weight):
        for i in range(weight):
            if self.nb_elmt<=1:
                if self.nb_elmt==1:
                    self.M = self.distribution.keys()[0]
                else:
                    self.M=0
                self.S = 0
            else:
                previous_M=self.M
                self.M = previous_M + ((value - previous_M) / self.nb_elmt)
                previous_S=self.S
                self.S = previous_S + (value - previous_M) * (value - self.M)

                
    def _update_std_dev_remove(self, value, weight):
        for i in range(weight):
            if self.nb_elmt<=1:
                if self.nb_elmt==1:
                    self.M = self.distribution.keys()[0]
                else:
                    self.M = 0
                self.S = 0
            else:
                previous_M=self.M
                self.M *= self.nb_elmt+1
                self.M = (self.M - value) / self.nb_elmt
                self.S -= (value - previous_M)*(value - self.M)
    
    def get_real_std_dev(self):
        """Calculate the standard deviation assuming the distribution is normal."""
        nb_element=self.get_nb_element()
        if nb_element<2:
            return 0
        sum2=0
        sumc=0
        mean=self.get_mean()
        set_of_value=self.distribution.keys()
        for val in set_of_value:
            sum2 = sum2 + ((val - mean)**2)*self.distribution.get(val)
            sumc = sumc + (val - mean)*self.distribution.get(val)
            
        variance = (sum2 - (sumc**2/nb_element))/(nb_element - 1)
        return math.sqrt(variance)
    
    def _get_std_dev_update(self):
        if self.nb_elmt>1:
            return math.sqrt(self.S / (self.nb_elmt - 1))
        else:
            return 0
    
    def get_percentiles(self,percentiles):
        """Calculate the percentiles of the distribution. 
        @param percentiles: a scalar or an array of scalar comprised in between 0 and 100.
        @return A value or a array value for the given percentiles."""
        set_of_value=self.distribution.keys()
        set_of_value.sort()
        if not percentiles.__class__ == [].__class__:
            percentiles=[percentiles]
        n_percentiles=[]
        sum=0
        for val in set_of_value:
            sum+=self.distribution.get(val)
        #print 'Sum=%s'%sum
        for percentile in percentiles:
            percentile_index=sum*(percentile/100.0)
            #print 'for %s percentile pos is %s'%(percentile,percentile_pos)
            index_sum=0
            if len(set_of_value)<1:
                n_percentiles.append(0)
            else:
                for val in set_of_value:
                    index_sum+=self.distribution.get(val)
                    if index_sum>=percentile_index:
                        #print 'val=%s\tsum2=%s: GOOD'%(val,sum2)
                        n_percentiles.append(val)
                        break
                    #print 'val=%s\tsum2=%s: NOT GOOD'%(val,sum2)
        if not percentiles.__class__ == [].__class__:
            return n_percentiles[0]
        else:
            return n_percentiles
    
    def get_greater_than(self, list_values):
        """Return the number of entries that are greater than the specified value"""
        nb_entries=self.get_nb_entry_with_generic_compare(list_values,lambda x,y: x>y)
        return nb_entries
    
    def get_greater_or_equal_than(self, list_values):
        """Return the number of entries that are greater or equal than the specified value"""
        nb_entries=self.get_nb_entry_with_generic_compare(list_values,lambda x,y: x>=y)
        return nb_entries
    
    def get_lower_than(self, list_values):
        """Return the number of entries that are lower than the specified value"""
        nb_entries=self.get_nb_entry_with_generic_compare(list_values,lambda x,y: x<y)
        return nb_entries
    
    def get_lower_or_equal_than(self, list_values):
        """Return the number of entries that are lower or equal than the specified value"""
        nb_entries=self.get_nb_entry_with_generic_compare(list_values,lambda x,y: x<=y)
        return nb_entries
    
    def get_nb_entry_with_generic_compare(self, list_values, compare):
        """Return a list of number of entries for which each element is the sum of all the entry that return true using the specified compare function."""
        set_of_value=self.distribution.keys()
        list_of_nb_entries=[0]*len(list_values)
        for tmp_value in set_of_value:
            for i,value in enumerate(list_values):
                if compare(tmp_value,value):
                    list_of_nb_entries[i]+=self.distribution.get(tmp_value)
        return list_of_nb_entries
    
    
    def get_new_distrib_with_generic_compare(self, compare, **xargs):
        """Return a distribution containing all the entries for which that return true using the specified compare function and specified value."""
        set_of_value=self.distribution.keys()
        new_distribution=Distribution_holder()
        for tmp_value in set_of_value:
            if compare(tmp_value,xargs):
                new_distribution.add_value(tmp_value, self.distribution.get(tmp_value))
        return new_distribution
    
    def get_cumsum(self):
        from numpy.ma.core import cumsum
        values=self.distribution.keys()
        values.sort()
        proportions=[]
        for value in values:
            proportions.append(self.distribution.get(value))
        cdf=cumsum(proportions)
        return cdf
            
    def get_cumsum_percent(self):
        from numpy.ma.core import cumsum
        values=self.distribution.keys()
        values.sort()
        proportions=[]
        for value in values:
            proportions.append(self.distribution.get(value))
        cdf=cumsum(proportions)
        total=sum(proportions)
        cdf_prop=[]
        for i in range(len(cdf)):
            cdf_prop.append(float(cdf[i])/total*100)
        return cdf_prop
    
    def plot(self, bins=None, output_file=None, xlimit=None, xlabel=None, ylabel=None, title=None, xlog=False,ylog=False):
        plot_distribution_holder(self, output_file=output_file, bins=bins,
                                 xlimit=xlimit, xlabel=xlabel, ylabel=ylabel,
                                 xlog=xlog,ylog=ylog,title=title)
    
    def print_dist(self,output_file=None,textgraph=None,sort_by_weight=False, reverse=False, nb_bin=None):
        return print_distribution_holder(self,output_file=output_file,textgraph=textgraph, sort_by_weight=sort_by_weight, reverse=reverse, nb_bin=nb_bin)
    
    def __str__(self):
        return print_distribution_holder(self,output_file=None)
    
    def add_distribution(self,dist):
        dict=dist.get_distribution()
        for value in dict.keys():
            weight=dict.get(value)
            self.add_value(value=value, weight=weight)
            
def load_distribution(open_stream, dist=None, string=False):
    if dist is None:
        dist=Distribution_holder()
    if string:
        func=str
    else:
        func=float
    for line in open_stream:
        sp_line=line.strip().split()
        if len(sp_line)>1:
            dist.add_value(func(sp_line[0]), int(sp_line[1]))
        if len(sp_line)>0:
            dist.add_value(func(sp_line[0]))
    return dist
            

def print_distribution_holder(holder, output_file=None, textgraph=None, sort_by_weight=False, reverse=False, nb_bin=None):
    if not sort_by_weight:
        if nb_bin is None:
            values,weights=holder.get_sorted_value_and_weight(reverse=reverse)
        else:
            values,weights=holder.get_binned_value(nb_bin=nb_bin)
    else:
        values,weights=holder.get_sorted_weight_and_value(reverse=reverse)
    
    out=[]
    if output_file:
        open_output=utils_logging.open_output_file(output_file, pipe=False)
        function=open_output.write
    else:
        function=out.append
    if textgraph:
        multiplier=200
        sum=0
        mark='|'
        maximum=max(weights)
        if maximum<multiplier:
            maximum=multiplier
        for i in range(len(values)):
            function('%s\t%s %s\n'%(values[i],(mark * int(float(weights[i])/maximum*multiplier)), weights[i] ))
    else:
        for i in range(len(values)):
            function('%s\t%s\n'%(values[i],weights[i]))
    if output_file:
        open_output.close()
        to_return = output_file
    else:
        to_return = ''.join(out)
    return to_return

def plot_distribution_holder(holder, output_file=None, bins=None, xlimit=None, xlabel=None, ylabel=None, title=None, xlog=False,ylog=False):
    try:
        import matplotlib.pyplot as plt
    except ImportError, e:
        logging.warning(str(e))
        logging.warning('You need to install matplotlib to use the graphical capabilities')
        return  
    if bins is None:
        if xlimit:
            new_holder = Distribution_holder()
            values, weights = holder.get_sorted_value_and_weight()
            for pos in range(len(values)):
                if values[pos]>xlimit[0] and values[pos]<xlimit[1]:
                    new_holder.add_value(values[pos], weights[pos])
        else:
            new_holder=holder
        percentile_90,=new_holder.get_percentiles([90])
        percentile_90,=new_holder.get_percentiles([100])
        percentile_0, = new_holder.get_percentiles([0])
        min_array=min(0, percentile_0)
        print min_array,float(percentile_90), float(percentile_90)/100
        bins=list(frange(min_array, float(percentile_90), increment = float(percentile_90)/100))
        if xlimit is None:
            xlimit=[bins[0],bins[-1]]
        
    ##Values for cdf
    cdf_prop=holder.get_cumsum_percent()
    
    #Values for histogram
    values, weights=holder.get_sorted_value_and_weight()
    fig = plt.figure()
    ax = fig.add_subplot(2,1,1)
    if title:
        plt.title(title)
    plt.hist(values, bins=bins, weights=weights)
    if xlog:
        ax.set_xscale('log')
    if ylog:
        ax.set_yscale('log')
    if not xlabel:
        xlabel='Value'
    plt.xlabel(xlabel)
    if not ylabel:
        tmp='number of %s'%(xlabel)
    else:
        tmp='number of %s'%(ylabel)
    plt.ylabel(tmp)
    if xlimit:
        plt.xlim(xlimit)
    
    ax2 = fig.add_subplot(2,1,2)
    plt.plot(values,cdf_prop)
    plt.xlabel(xlabel)
    if not ylabel:
        tmp='proportion of %s (%%)'%(xlabel)
    else:
        tmp='proportion of %s (%%)'%(ylabel)
    plt.ylabel(tmp)
    if xlimit:
        plt.xlim(xlimit)
    plt.ylim([0,100])
    if output_file:
        plt.savefig(output_file)
    else:
        plt.show()
        

class Binned_Distribution_Holder(object):
    
    def __init__(self, bin_size):
        self.all_bins={}
        self.bin_size
    
    def add_pair(self, key, value):
        dist=self.all_bins.get(int(key/self.bin_size))
        if dist is None:
            dist = Distribution_holder()
            self.all_bins[int(key/self.bin_size)]=dist
        dist.add_value(value)
    
    def get_binned_distribution(self):
        return self.all_bins
    
    
    
class Histogram(object):
    '''
    classdocs
    '''


    def __init__(self, bin_size):
        '''
        Constructor
        '''
        self.bin_size=bin_size
        

def bin_value_from_array(array, bin_size, max_value=None):
    all_bins=[]
    if max_value:
        if type(bin_size) == float:
            range=frange
        for i in range(0,max_value,bin_size):
            all_bins.append(0)
    for value in array:
        bin = value/bin_size
        if len(all_bins)<=bin:
            to_add=bin-len(all_bins)+1
            for i in range(to_add):
                all_bins.append(0)
        all_bins[bin]+=1
    return all_bins

def bin_coordinates(input_file, output_file, bin_size):
    open_file=utils_logging.open_input_file(input_file)
    open_output=utils_logging.open_output_file(output_file)
    all_coordinates_per_chr={}
    for line in open_file:
        sp_line=line.split()
        all_coordinates=all_coordinates_per_chr.get(sp_line[0])
        if all_coordinates is None:
            all_coordinates=[]
            all_coordinates_per_chr[sp_line[0]]=all_coordinates
        all_coordinates.append(int(sp_line[1]))
    
    for chr in all_coordinates_per_chr.keys():
        all_coordinates=all_coordinates_per_chr.get(chr)
        all_bins=bin_value_from_array(all_coordinates, bin_size)
        for bin,value in enumerate(all_bins):
            open_output.write('%s\t%s\t%s\n'%(chr,bin*bin_size,value))
    open_output.close()


def bin_coordinates_through_genome(input_file, output_file, genome_file, bin_size):
    open_file=utils_logging.open_input_file(input_file)
    open_output=utils_logging.open_output_file(output_file)
    all_coordinates_per_chr={}
    genome_loader=GenomeLoader(genome_file)
    previous_bin=0
    all_chr=[]
    for line in open_file:
        sp_line=line.split()
        all_coordinates=all_coordinates_per_chr.get(sp_line[0])
        if all_coordinates is None:
            all_chr.append(sp_line[0])
            all_coordinates=[]
            all_coordinates_per_chr[sp_line[0]]=all_coordinates
        all_coordinates.append(int(sp_line[1]))
    all_chr.sort()
    for chr in all_chr:
        header, sequence =genome_loader.get_chr(chr)
        chr=header.strip()
        chr_len=len(sequence)
        
        all_coordinates=all_coordinates_per_chr.get(chr)
        all_bins=bin_value_from_array(all_coordinates, bin_size, chr_len)
        for bin,value in enumerate(all_bins):
            open_output.write('%s\t%s\t%s\t%s\n'%(chr, bin*bin_size, (bin*bin_size)+previous_bin, value))
        previous_bin+=len(all_bins)*bin_size
    open_output.close()    
    
if __name__=="1__main__":
    input_file='/home/tcezard/projects/2009127_carol_wanqui/provided_snps_list/3D7HB3_sorted_solexa_data_filtered_40_coverage.csv'
    output_file='/home/tcezard/projects/2009127_carol_wanqui/provided_snps_list/3D7HB3_sorted_solexa_data_filtered_40_coverage_bin.csv'
    genome_file='/home/tcezard/genomes/plasmodium_falciparum/PlasmoDB-6.3/PfalciparumGenomic_PlasmoDB-6.3.fasta'
    bin_size=10000
    bin_coordinates_through_genome(input_file, output_file, genome_file, bin_size)


   
if __name__=="1__main__":
    import random
    dist=Distribution_holder(update=False)
    mean=1000
    sigma=10
    all_values=[]
    difference_add=[]
    difference_remove=[]
    sample_size=1000
    for i in range(sample_size):
        value=random.normalvariate(mean, sigma)
        dist.add_value(value)
        all_values.append(value)
        std_dev=dist.get_std_dev()
        real_std_dev=dist.get_real_std_dev()
        print abs(real_std_dev-std_dev)
        difference_add.append(abs(real_std_dev-std_dev))
        
    for value in all_values:
        dist.remove_value(value)
        std_dev=dist.get_std_dev()
        real_std_dev=dist.get_real_std_dev()
        print abs(real_std_dev-std_dev)
        difference_remove.append(abs(real_std_dev-std_dev))
    difference_remove=difference_remove[::-1]
    x=range(1,sample_size+1)
    print len(x), len(difference_add), len(difference_remove)
    import matplotlib.pyplot as plt
    plt.plot(x,difference_add, '-b', linewidth=2, label='adding number')
    plt.plot(x,difference_remove,'-r', label='removing number')
    plt.xlabel('nb element in the distribution')
    plt.ylabel('abs(real std dev - updated std dev)')
    plt.legend()
    
    plt.show()

if __name__=="__main__":
    import random, math
    dist=Distribution_holder(update=False)
    dist2=Distribution_holder(update=True)
    mean=450
    sigma=10
    sample_size=1000
    test_dist=[dist,dist2]
    for d in test_dist:
        for i in range(sample_size):
            value=random.normalvariate(mean, sigma)
            d.add_value(value)
            number_element=d.get_nb_element()
        try:
            assert number_element==sample_size
        except AssertionError:
            print 'number element in distribution %s is different from number of element entered %s'%(number_element,sample_size)
        try:
            assert round(d.get_mean()) == mean
        except AssertionError:
            print 'mean in distribution %s is different from the mean use to draw random number %s'%(d.get_mean(),mean)
        try:
            assert round(d.get_std_dev()) == sigma
        except AssertionError:
            print 'Standard deviation in distribution %s is different from the standard deviation use to draw random number %s'%(d.get_std_dev(),sigma)
        
        
