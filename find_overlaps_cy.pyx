import pandas as pds
import numpy as np
import math


def create_overlap_table(list_arr_df, min_ovl_len, path1):
    arr_1,start_df = list_arr_df
    
    overlap_list = []
    
    #find overlaps in sorted list of reads and append to overlap_list
    for index,list1 in enumerate(arr_1):
        for list2 in arr_1[index+1:]:
            if list2[1] < list1[2]:
                overlap_list.append(list1[0:1] + list2)
            else:
                break
    
    cols_to_exclude = ("start_pos end_minus_ovl").split()
    all_cols = ("query_id1 query_id2").split() + cols_to_exclude
    
    df0 = pds.DataFrame.from_records(overlap_list, exclude=cols_to_exclude, columns=all_cols) #convert overlap_list to pandas df
    
    df0 = df0.merge(start_df, left_on=all_cols[0], right_on='query_id', how='inner').drop('query_id', axis=1)

    df0 = df0.merge(start_df, left_on=all_cols[1], right_on='query_id', how='inner', suffixes=('1','2')).drop('query_id', axis=1)
    
    cols_to_remove = ("start_pos1 start_pos2 end_pos1 end_pos2").split()
    final_cols = ("query_id1 query_id2 read_length1 read_length2 overlap_length strand_orient1 strand_orient2").split()
    
    df0['overlap_length'] = df0[cols_to_remove[2:]].min(axis=1) - df0[cols_to_remove[1]]
    
    df0 = df0.drop(cols_to_remove, axis=1)[final_cols].sort_values(final_cols[0:2])
    
    df0.to_csv(path1, header=False, sep='\t', index=False)







def check_data(dataframe1, overlap_length):
    #sanity checks
    if overlap_length < 1:
        raise ValueError('overlap_length must be >= 1') #raise error if overlap_length is less than 1
    
    if type(overlap_length) != int:
        raise TypeError('overlap_length must be int') #raise error if overlap_length is not an int 
    
    if dataframe1.shape[1] != 6:
        raise ValueError('BLAST output must have 6 columns: query_id, start_pos, end_pos, e_value, bit_score, and percent_identity') #raise error if there are not 6 columns in BLAST data
        
    col_names = ['query_id','start_pos','end_pos','e_value','bit_score','p_ident'] 
    dataframe1.columns = col_names  
    
    if any(dataframe1['start_pos'] <= 0) or any(dataframe1['end_pos'] <= 0) or any(dataframe1['e_value'] < 0) or any(dataframe1['p_ident'] < 0) or any(dataframe1['p_ident'] > 100):
        raise ValueError('invalid values')

    if dataframe1['query_id'].dtype == 'object':
        try:
            dataframe1['query_id'] = dataframe1['query_id'].astype('int64') 
        except ValueError:
            dataframe1['query_id'] = dataframe1['query_id'].map(lambda x: x[x.index('.') + 1:]).astype('int64') #if query_id is not a string with all numeric values, remove non-numeric 
                                                                                                                #characters from and query_id and convert to int
    
    if any(dataframe1.duplicated('query_id',False)):
        sub_dataframe1 = dataframe1.loc[dataframe1.duplicated('query_id', False), col_names]
        
        dataframe1_no_dup = dataframe1.drop_duplicates('query_id').set_index('query_id')
        
        sub_dataframe1_no_dup = sub_dataframe1.sort_values(['query_id','e_value','bit_score'], 
                                    ascending=[True,True,False]).drop_duplicates('query_id',
                                                                                'first').set_index('query_id') #sort each duplicated query_id by e_value and bit_score,then drop
                                                                                                               #all duplicates except first duplication in sorted dataframe
        no_dup_q_id = sub_dataframe1_no_dup.index
        dataframe1_no_dup.loc[no_dup_q_id] = sub_dataframe1_no_dup
    
    dataframe1_no_dup = dataframe1_no_dup.drop(('e_value bit_score p_ident').split(), axis=1)
    dataframe1_no_dup['read_length'] = (dataframe1_no_dup['end_pos'] - dataframe1_no_dup['start_pos'])
    dataframe1_no_dup['strand_orient'] = np.where(dataframe1_no_dup['read_length']<0, 'R', 'F') #specify strand type, reverse or forward, based on whether read_length is positive (F) or negative (R)
    dataframe1_no_dup['read_length'] = dataframe1_no_dup['read_length'].abs() #make sure that read_length is positive
    dataframe1_no_dup['read_length'] = dataframe1_no_dup['read_length'] + 1
    
    np_arr = dataframe1_no_dup.loc[:,['start_pos','end_pos']].values
    np_arr.sort(axis=1)
    dataframe1_no_dup.loc[:,['start_pos','end_pos']] = pds.DataFrame(np_arr, dataframe1_no_dup.index, ['start_pos','end_pos']) #make start_pos < end_pos
    
    dataframe1_no_dup = dataframe1_no_dup.sort_values(['start_pos', 'end_pos'], ascending=[True,True]).reset_index() #sort by start_pos, then by end_pos, both in ascending order
    
    dataframe1_no_dup = dataframe1_no_dup.loc[dataframe1_no_dup['read_length'] >= overlap_length] #filter out reads with read_length less than overlap_length


    dataframe1_no_dup['end_minus_ovl'] = (dataframe1_no_dup['end_pos'] - overlap_length) + 1
    cols_to_keep = ("query_id start_pos end_minus_ovl").split()
    arr = dataframe1_no_dup[cols_to_keep].to_numpy(dtype='i8',copy=True).tolist()
    
    dataframe1_no_dup.drop('end_minus_ovl', inplace=True, axis=1)
    
    list0 = [arr, dataframe1_no_dup]

    return list0
    
    
    
    
    
def find_overlaps(filepath1, savepath1, overlap_length1=50):
    df = pds.read_csv(filepath1, header=None, sep='\t')
    
    list_0 = check_data(df, overlap_length1)
    
    create_overlap_table(list_0, overlap_length1, savepath1)
    
	
	
	
	

#find_overlaps(filepath1='/user/jzola/projects/mepiskoz/identify overlaps/ERR2906227_blast/ERR2906227_g10.blast', overlap_length1 = 10, 
#savepath1='/user/jzola/projects/mepiskoz/identify overlaps/overlaps/overlaps_len10_g10')





