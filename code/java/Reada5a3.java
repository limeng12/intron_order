import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.bed.FullBEDFeature.Exon;
import htsjdk.tribble.readers.LineIterator;

public class Reada5a3 {

	public static void main(String[] args) throws Exception {
		//read5a3("");
		
		Cigar c=new Cigar();
		//c.add(new CigarElement(16,CigarOperator.S) );

		//c.add(new CigarElement(10,CigarOperator.D) );
		//c.add(new CigarElement(12,CigarOperator.I) );

		//c.add(new CigarElement(13,CigarOperator.N) );
		//c.add(new CigarElement(16,CigarOperator.S) );
		c.add(new CigarElement(257,CigarOperator.M) );
		c.add(new CigarElement(57,CigarOperator.N) );
		c.add(new CigarElement(136,CigarOperator.M) );
		c.add(new CigarElement(168,CigarOperator.N) );
		c.add(new CigarElement(180,CigarOperator.M) );

		
		
		int [] aa=get_splice_junction_pos(c,180272);

		
		//int [] aa=get_splice_junction_pos(c,1);
		
		for( int i =0;i< aa.length;i++) {
			
			System.out.println(aa[i]);
		}
		
		
	}	
	
	public static HashMap<String, Integer> read5a3(){
		return read5a3("");

	}
	
	public static HashMap<String, Integer> read5a3_no5a3(){
		//return read5a3("");
		HashMap<String, Integer> e_map = new HashMap<String, Integer>();
		
		return e_map;
	}
	
	
	public static HashSet<String> readTranscripts(String t_path) {
		HashSet<String> m_trans = new HashSet<String>();
		
		try {
		BufferedReader buf;

		buf = new BufferedReader(new FileReader(t_path) );
		
		String lineJustFetched;
		String[] wordsArray;

		while (true) {
			lineJustFetched = buf.readLine();
			if (lineJustFetched == null) {
				break;
			} else {
				//wordsArray = lineJustFetched.split("\t");
				m_trans.add(lineJustFetched);
				//e_map.put((wordsArray[0]), Integer.parseInt(wordsArray[1]));
				//System.out.println(wordsArray[0]);
				
			}
		}
		
		buf.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return m_trans;
		
	}
	
	public static HashMap<String, Integer> read5a3(String t_path) {
		HashMap<String, Integer> e_map = new HashMap<String, Integer>();

		String p=Reada5a3.class.getResource("/resources/a3a5_bias.tsv").getPath();
		
		try {
			BufferedReader buf;
			
			String path="";
			
			if(t_path=="") {
				path=p;
			}else {
				path=t_path;
			}
			
			
			System.out.println("a3a5_bias_file: "+path);
			
			buf = new BufferedReader(new FileReader(path) );

			
			String lineJustFetched = null;
			String[] wordsArray;

			while (true) {
				lineJustFetched = buf.readLine();
				if (lineJustFetched == null) {
					break;
				} else {
					wordsArray = lineJustFetched.split("\t");
					e_map.put((wordsArray[0]), Integer.parseInt(wordsArray[1]));
					//System.out.println(wordsArray[0]);
					
				}
			}
			
			buf.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
		return e_map;

	}
	
	public static int convert_cigar_to_len(CigarElement ele) {
		
		int cover_len=0;
		
        switch (ele.getOperator()) {
	        case M:
	        case D:
	        case N:
	        case EQ:
	        case X:
	        	cover_len = ele.getLength();
	            break;
	        default: break;
	    }
		
		
		return cover_len;
	}
	
	public static HashSet<String> get_all_intron(File t_bed_file) throws IOException{
		
		AbstractFeatureReader<BEDFeature, LineIterator> br = AbstractFeatureReader.getFeatureReader(t_bed_file.getPath(),
				new BEDCodec(), false);
		
		Iterable<BEDFeature> iter = br.iterator();

		HashSet<String> all_introns = new HashSet<String>();
		
		
		for (BEDFeature feat : iter) {
			
			int count_of_exon = feat.getExons().size();

			for (int i = 1; i < (count_of_exon ); i++) {

				Exon e=feat.getExons().get(i);

				String chr =feat.getChr();

				int start = e.getCdStart();

				int left_int_coor=feat.getExons().get(i-1).getCdEnd()+1;

				String left_int = chr+":"+left_int_coor+"-"+(start-1);
				
				all_introns.add(left_int);
				
			}
			
		}
		
		br.close();
		
		return all_introns;
	}
	
	
	public static boolean all_junctions_in_intron(HashSet<String> t_all_introns,HashSet<String> t_all_introns_trans,String t_chr, int[] t_star_ss) {
		for ( int j = 0; j < (Math.floor(t_star_ss.length / 2) ); j++) {
			
			//String one_region = recode.getContig() + ":" + (att[1 + i * 2] + 1) + "-"
			//		+ (att[2 + i * 2] - 1);
			//String one_region = recode.getContig() + ":" + (att[ j * 2] ) + "-"+ (att[ j * 2+1] );
			
			String one_region=t_chr+":"+(t_star_ss[ j * 2]) +"-" + (t_star_ss[ j * 2+1] ) ;
			
			/*
			 * if junctions in a read not in the transcript's intron, but in this gene, then this read is not compatible with 
			 * the transcript, more usefull of this function if longer read
			 */
			if( (!t_all_introns_trans.contains(one_region) ) & (t_all_introns.contains(one_region))   ) {
				return false;
			}
			
		}
		
		
		return true;
		
	}
	
	
	public static int[] get_splice_junction_pos(Cigar t_cigar,int start) {
		
		ArrayList<Integer> jc_pos=new ArrayList<Integer>();
		
		int current_len=start;
		
		for(CigarElement tt_c : t_cigar.getCigarElements()) {
			CigarOperator t_ope=tt_c.getOperator();

			int cov_len=convert_cigar_to_len(tt_c );
			
			if(t_ope.equals(CigarOperator.N) & tt_c.getLength()>10 ) {
				jc_pos.add(current_len);
				jc_pos.add(current_len+ tt_c.getLength()-1);
				//jc_pos.add(current_len);
				//jc_pos.add(current_len+ tt_c.getLength()-1);
			}
			
			current_len=current_len+cov_len;
		}
		
		
		int[] jc_pos_arr = new int[jc_pos.size()];
		
		for(int i=0;i<jc_pos_arr.length;i++) {
			jc_pos_arr[i]=jc_pos.get(i);
		}
		
				
		return jc_pos_arr;
		
		
	}
	
	
	
}


