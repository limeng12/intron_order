import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.bed.FullBEDFeature.Exon;
import htsjdk.tribble.readers.LineIterator;

public class ISOLarge {

	public static void main(String[] args) throws Exception {
		HashMap<String, Integer> e_map;
		
		
		/*
		//String in_bed_path = "/Users/mengli/Documents/projects/iso/anno/hg19_gencode_from_ucsc_nothick_nocds_test.bed";
		String in_bed_path = "/Users/mengli/Documents/projects/iso/anno/hg19_gencode_from_ucsc_nothick_nocds.bed";

		String bam_path = "/Users/mengli/Downloads/ENCFF804MYH.bam";
		
		String out_path = "/Users/mengli/Documents/projects/iso/result/iso_test.tsv";

		int intron_flank_threshold=20;
		
		boolean unique_int=true;
		
		boolean consider_large_int=true;
		
		boolean consider_exon_in_intron=true;
		*/
		
		
		
		/*
		String in_bed_path = "/Users/mengli/Documents/projects/iso/zebrafish/GRCz11_ensembl_nothick.bed";

		String bam_path = "/Users/mengli/Downloads/";
		
		String out_path = "/Users/mengli/Documents/projects/iso/result/iso_test.tsv";

		int intron_flank_threshold=20;
		
		boolean unique_int=true;
		
		boolean consider_large_int=true;
		
		e_map=Reada5a3.read5a3("/Users/mengli/Documents/projects/iso/zebrafish/a3a5_bias_all_zebra.tsv");
		 */
		
		String in_bed_path=args[0];
		
		String bam_path=args[1];
		
		String out_path = args[2];
		
		/*
		HashSet<String> target_trans=new HashSet<String>();
		
		if(args.length>3) {
			
			target_trans=Reada5a3.readTranscripts(args[3]);
		
		}
		*/
		
		int intron_flank_threshold=20;
		
		if(args.length>3) {
			intron_flank_threshold=Integer.parseInt(args[3]);
		}
		
		boolean unique_int=false;
		if(args.length>4) {
			unique_int=Boolean.parseBoolean(args[4]);
		}
		
		boolean consider_large_int=true;
		if(args.length>5) {
			consider_large_int=Boolean.parseBoolean(args[5]);
		}
		
		
		//one isoform's exon maybe in another isoform's intron region
		boolean consider_exon_in_intron=true;
		
		if(args.length>6) {
			e_map=Reada5a3.read5a3(args[6]);
		}else {
			//e_map=Reada5a3.read5a3();
			e_map=Reada5a3.read5a3_no5a3();
		}	

		if(args.length>7) {
			consider_exon_in_intron=false;
		}
		
		System.out.println("consider whole transcript: "+consider_large_int);
		
		System.out.println("unique intron: "+unique_int);
		
		
		
		File bed_file = new File(in_bed_path);
		
		// BufferedReader br = new BufferedReader(new FileReader(file));
		
		AbstractFeatureReader<BEDFeature, LineIterator> br = AbstractFeatureReader.getFeatureReader(bed_file.getPath(),
				new BEDCodec(), false);
		
		
		File bamFile = new File(bam_path);
		//File baiFileIndex = new File(bamFile.getPath() + ".bai");
		
		// private BAMFileReader bamFileReaderBAI = new BAMFileReader(bamFile,
		// baiFileIndex, true, false,
		// ValidationStringency.DEFAULT_STRINGENCY,
		// DefaultSAMRecordFactory.getInstance());

		SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
		// File output_one= new File("output");
		BufferedWriter output_one = new BufferedWriter(new FileWriter(out_path));

		HashMap<String, Integer> interest_iso = new HashMap<String, Integer>();

		//String st;
		Iterable<BEDFeature> iter = br.iterator();

		HashSet<String> two_int_coord = new HashSet<String>();
		
		HashSet<String> all_introns=Reada5a3.get_all_intron(bed_file);
		
		
		for (BEDFeature feat : iter) {

			// System.out.println(st);
			int count_of_exon = feat.getExons().size();
			if (count_of_exon < 3) {
				continue;
			}
			
			
			
			if(!feat.getChr().startsWith("chr") ) {
				continue;
			}
			
			
			if(feat.getChr().contains("gl") | feat.getChr().contains("Un_") ) {
				continue;
			}
			
			String strand=feat.getStrand().toString();
			String trans_name=feat.getName();
			
			//String[] array = trans_name.split("\\.",-1); 
			
			//trans_name=array[0];
			
			//if(!target_trans.contains(trans_name) ) {
			//	continue;
			//}
			
			
			int trans_start=feat.getStart();
			
			int trans_end=feat.getEnd();
			
			

			/*
			 * record all intron coordiantes, if juntioin containing read doesn't compatible with 
			 * the intron, then they are not in this transcript. 
			 */
			HashSet<String> all_intron_trans=new HashSet<String>();
			
			for (int i = 1; i < (count_of_exon ); i++) {

				Exon e=feat.getExons().get(i);

				String chr =feat.getChr();

				int start = e.getCdStart();

				int left_int_coor=feat.getExons().get(i-1).getCdEnd()+1;

				String left_int = chr+":"+left_int_coor+"-"+(start-1);
				
				all_intron_trans.add(left_int);
				
			}
			
			
			
			for (int i = 1; i < (count_of_exon ); i++) {

				Exon e=feat.getExons().get(i);
				
				String chr =feat.getChr();
				int start = e.getCdStart();
				
				//int end = e.getCdEnd();
				
				
				
				int left_int_coor=feat.getExons().get(i-1).getCdEnd()+1;

				String left_int = chr+":"+left_int_coor+"-"+(start-1);
				
				
				HashSet<String> cover_left_names = new HashSet<String>();
				
				//HashSet<String> cover_right_names = new HashSet<String>();
				
				if(unique_int & two_int_coord.contains(chr+":"+left_int) ) {
					continue;
				}
				
				two_int_coord.add(chr+":"+left_int);
				
				
				
				int right_exon_start=start;
				
				int left_exon_end=left_int_coor-1;
				
				boolean exon_ri_a3_a5_likely_bias=false;
				
				//correct for a5ss and a3ss, to avoid bias of exon consider as introns
				if(e_map.containsKey(chr+":"+start)) {
					int t_bound_start=e_map.get(chr+":"+start);
					if(t_bound_start>left_int_coor) {
						right_exon_start=t_bound_start;
					}else {
						exon_ri_a3_a5_likely_bias=true;
					}
					
				}
				
				//correct for a5ss and a3ss, to avoid bias of exon consider as introns
				if(e_map.containsKey(chr+":"+(left_exon_end) )) {
					int t_left_int_coor=e_map.get(chr+":"+(left_exon_end));
					if(right_exon_start>t_left_int_coor) {
						left_exon_end=t_left_int_coor;
					}else {
						exon_ri_a3_a5_likely_bias=true;

					}
					
				}
				
				
				//CloseableIterator<SAMRecord> r = reader.queryContained(chr, start, end);
				
				//System.out.println(chr+":"+left_int_coor+"-"+(start-1));
				CloseableIterator<SAMRecord> r = null;
				//System.out.println(chr);System.out.println(left_int_coor);System.out.println(start);

				try {
					r=reader.query(chr, left_int_coor, start-1,false);
				}catch(Exception e1) {
					System.out.println(e1.getMessage());
					e1.printStackTrace();
					System.out.println(chr);System.out.println(left_int_coor);System.out.println(start);

					continue;
				}
				//reader.queryContained(arg0, arg1, arg2)
				while (r.hasNext()) {
					
					SAMRecord recode = r.next();
					// String att=recode.getStringAttribute("jI");
					// boolean is_neg = recode.getReadNegativeStrandFlag();
					
					//bwa, mapping quality >0 is unique map
					if( (recode.getMappingQuality()<1)  ) {
						continue;
					}
					
					//star mapping, 255 is unique mapping
					//if( (recode.getMappingQuality()<255) & recode.hasAttribute("jI") ) {
					//	continue;
					//}
					
					if(recode.getReadLength()<30 ) {
						continue;
					}
					
					//STAR has NH tag, >1 not unique, =1 unique mapping
					if(recode.hasAttribute("NH") ) {
						if(recode.getIntegerAttribute("NH")>1) {
							continue;
						}
						
					}
					
					
					//boolean if_first=recode.getFirstOfPairFlag();
					//String name = recode.getReadName()+if_first;
					String name = recode.getReadName();

					// String pair_name=recode.getPairedReadName();

					
					//boolean cover_left = (recode.getStart() < (start - 10)) & (recode.getEnd() > (start + 10));
					//boolean cover_right = (recode.getStart() < (end - 10)) & (recode.getEnd() > (end + 10));

					
					
					boolean cover_left=false;
					//boolean cover_right=false;
					
					boolean has_N = recode.getCigar().containsOperator(CigarOperator.N);
					
					int anchor_region_len=100;
					
					if(anchor_region_len> (right_exon_start-left_exon_end-1)) {
						anchor_region_len=right_exon_start-left_exon_end-1;
					}
					
					if(intron_flank_threshold> (right_exon_start-left_exon_end-1)) {
						intron_flank_threshold=right_exon_start-left_exon_end-1-1;
					}
					
					
					if(has_N) {
						//int[] att = recode.getSignedIntArrayAttribute("jI");
						
						int[] att =null;
						
						if(recode.hasAttribute("jI")) {
							att=recode.getSignedIntArrayAttribute("jI");
						}else {
							att=Reada5a3.get_splice_junction_pos(recode.getCigar(),recode.getStart() );
						}
						
						
						if(att.length==0) {
							continue;
						}
						
				        //int att_min = Arrays.stream(att).min().getAsInt();
				        //int att_max = Arrays.stream(att).max().getAsInt();
						
						/*
						 * has junction, and junction are compatible with transcript
						 */
				        if(att.length>0) {
							boolean is_all_junctions_in_intron=Reada5a3.all_junctions_in_intron(all_introns,all_intron_trans, chr, att);
					        
							if(!is_all_junctions_in_intron) {
								continue;
							}
				        }
						
				        boolean no_left_sj=true;
				        //boolean no_right_sj=true;
				        
				        //no_left_sj=(att_min > start-1) | (att_max < left_int_coor);
				        
				        //no junction site in the intron
				        for(int tt=0;tt<att.length;tt++) {
				        	if( (att[tt]<=start-1) & (att[tt]>=left_exon_end+1 ) ) {
				        		no_left_sj=false;
				        	}
				        	
				        }
				        
				        
				        //no overlap junc span the intron,
				        //if no junc span the intron, either, the junction are in right or left of the intron
				        // if junc in both the right and left intron, then the read must cover the whole intron
				        if(no_left_sj==true) {
							for ( int j = 0; j < (Math.floor(att.length / 2) ); j++) {
	
								//String one_region = recode.getContig() + ":" + (att[1 + i * 2] + 1) + "-"
								//		+ (att[2 + i * 2] - 1);
								//String one_region = recode.getContig() + ":" + (att[ j * 2] ) + "-"+ (att[ j * 2+1] );
								if((att[ j * 2]<=start-1) &  (att[ j * 2+1]>=left_int_coor ) ) {
					        		no_left_sj=false;
					        		continue;
								}
							
							}
				        }
				        
						
						
						if(consider_exon_in_intron) {
				        /*
				         *  force the read overlap with intron at least intron_flank_threshold bp. 
				         *  read must around with intron-exon junction ( -80bp->0) to ensure not spliced. 
				         *  and the minimum bases in this region is intron_flank_threshold
				         *  
				         */
							
						cover_left = (no_left_sj) & (
								( (recode.getStart() <= (right_exon_start - intron_flank_threshold)) &
								(recode.getEnd() >= (right_exon_start - anchor_region_len+intron_flank_threshold) )    ) |
								(   (recode.getStart() <=(left_exon_end + anchor_region_len-intron_flank_threshold) )  &
								(recode.getEnd() >= (left_exon_end + intron_flank_threshold)) )
								);
						}else {
							cover_left = (no_left_sj) & (
									( (recode.getStart() <= (right_exon_start - intron_flank_threshold))    ) &
									(   
									(recode.getEnd() >= (left_exon_end + intron_flank_threshold)) )
									);
						}
						
						
				        /*
						cover_left = (end<att_min) & (
								( (recode.getStart() < (bound_start - intron_flank_threshold)) &
								(recode.getEnd() > (bound_start - 80) )    ) |
								(   (recode.getStart() < (left_int_coor + 80) )  &
								(recode.getEnd() > (left_int_coor + intron_flank_threshold)) )
								);
						*/					
				 
					}else {
						if(consider_exon_in_intron) {

						cover_left = ( (recode.getStart() <= (right_exon_start - intron_flank_threshold)) &
								(recode.getEnd() >= (right_exon_start - anchor_region_len+intron_flank_threshold) )    ) |
								(   (recode.getStart() <= (left_exon_end + anchor_region_len-intron_flank_threshold) )  &
								(recode.getEnd() >= (left_exon_end + intron_flank_threshold)) );
						}else {
							
							cover_left = ( (recode.getStart() <= (right_exon_start - intron_flank_threshold))  ) &
									(  
									(recode.getEnd() >= (left_exon_end + intron_flank_threshold)) );
						}
				
						
					}
					
					if (cover_left)
						cover_left_names.add(name);
					/*
					if (cover_right)
						cover_right_names.add(name);
					*/
					
				}
				
				r.close();
				
				CloseableIterator<SAMRecord> r_sec = null;
				
				if(cover_left_names.size()==0) {
					continue;
				}
				
				if(consider_large_int) {
					r_sec = reader.query(chr, trans_start, trans_end,false);
				}else {
					r_sec = reader.query(chr, trans_start, trans_end,false);
				}
				
				
				while (r_sec.hasNext()) {

					SAMRecord recode = r_sec.next();
					
					if( (recode.getMappingQuality()<1)  ) {
						continue;
					}
					
					//star mapping, 255 is unique mapping
					//if( (recode.getMappingQuality()<255) & recode.hasAttribute("jI") ) {

					//	continue;
					//}
					
					if(recode.getReadLength()<30 ) {

						continue;
					}
					
					//STAR has NH tag, >1 not unique, =1 unique mapping
					if(recode.hasAttribute("NH") ) {
						if(recode.getIntegerAttribute("NH")>1) {
							continue;
						}
						
					}

					/*
					if (!recode.getProperPairFlag() ) {
						continue;
					}
					 */
					
					
					int[] att =null;
					
					if(recode.hasAttribute("jI")) {
						att=recode.getSignedIntArrayAttribute("jI");
					}else {
						att=Reada5a3.get_splice_junction_pos(recode.getCigar(),recode.getStart() );
					}
					
					boolean has_jc = att.length > 0;

					if (!has_jc) {

						continue;
					}
					
					
					boolean is_all_junctions_in_intron=Reada5a3.all_junctions_in_intron(all_introns,all_intron_trans, chr, att);
			        
					if(!is_all_junctions_in_intron) {
						continue;
					}
					
					
					
					String low_int = "";

					//boolean if_first=recode.getFirstOfPairFlag();
					
					//String paired_name=recode.getReadName()+!if_first;
					String paired_name=recode.getReadName();

					//String paired_name=recode.getPairedReadName();
					
					//String paired_name=name;
					
					if (cover_left_names.contains(paired_name) ) {
						low_int = left_int;
					} else {
						continue;
					}
					
					if (att.length < 2) {

						continue;
					}

					// for(int i=0;i<((att.length/2)-1);i++) {
					for ( int j = 0; j < (Math.floor(att.length / 2) ); j++) {

						//String one_region = recode.getContig() + ":" + (att[1 + i * 2] + 1) + "-"
						//		+ (att[2 + i * 2] - 1);
						String one_region = recode.getContig() + ":" + (att[ j * 2] ) + "-"+ (att[ j * 2+1] );
						
						one_region = trans_name+"\t"+low_int + "\t" + one_region+"\t"+strand+"\t"+exon_ri_a3_a5_likely_bias;

						if (!interest_iso.containsKey(one_region)) {
							interest_iso.put(one_region, 1);
							continue;
						}

						int one_rcount = interest_iso.get(one_region);
						interest_iso.put(one_region, one_rcount + 1);

						//output_one.write(m.getKey() + "\t" + m.getValue());
						//output_one.write(output_one+"\t"+(one_rcount+1) );
						//output_one.newLine();

					}

					/*
					 * if(att.length>4) { String
					 * one_region=recode.getContig()+":"+(att[1]+1)+"-"+(att[2]-1);
					 * all_exon.add(one_region);
					 * 
					 * 
					 * one_region=recode.getContig()+":"+(att[3]+1)+"-"+(att[4]-1);
					 * all_exon.add(one_region); //System.out.println( one_region ); continue; }
					 * 
					 * 
					 * if(att.length>2) { String
					 * one_region=recode.getContig()+":"+(att[1]+1)+"-"+(att[2]-1);
					 * all_exon.add(one_region); //System.out.println( one_region ); continue; }
					 */
					
				}
				r_sec.close();
				
			}
			
			
			for (Entry<String, Integer> m : interest_iso.entrySet()) {
				output_one.write(m.getKey() + "\t" + m.getValue());
				output_one.newLine();
			}
			
			interest_iso.clear();
			
		}

		/*
		for (Entry<String, Integer> m : interest_iso.entrySet()) {
			output_one.write(m.getKey() + "\t" + m.getValue());
			output_one.newLine();

		}
		*/
		
		for (Entry<String, Integer> m : interest_iso.entrySet()) {
			output_one.write(m.getKey() + "\t" + m.getValue());
			output_one.newLine();
		}
		
		interest_iso.clear();
		
		output_one.close();

		/*
		 * it = all_exon.iterator(); // why capital "M"? while(it.hasNext()) {
		 * output_one.write(it.next()); output_one.newLine(); } output_one.close();
		 */

	}
	
}
