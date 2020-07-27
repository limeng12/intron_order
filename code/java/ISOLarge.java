import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.bed.FullBEDFeature.Exon;
import htsjdk.tribble.readers.LineIterator;

import org.apache.commons.cli.*;

public class ISOLarge {

	public static HashSet<String> preprocess_name(SamReader t_bam_reader) throws IOException {

		
		SAMFileHeader t_bam_header = t_bam_reader.getFileHeader();
		
		SAMSequenceDictionary bam_header_dic = t_bam_header.getSequenceDictionary();
		
		HashSet<String> seq_names=new HashSet<String>();
		
		for( SAMSequenceRecord r:bam_header_dic.getSequences() ) {
			
			String seq_name_one=r.getSequenceName();
			seq_names.add(seq_name_one);
			
			bam_header_dic.addSequenceAlias(seq_name_one, "chr"+seq_name_one);
			//System.out.println(seq_name_one);
			
		}
		
				
		return seq_names;
		
	}
	
	public static void test_iso() throws Exception {

		
		//String in_bed_path = "/Users/mengli/Documents/projects/iso/anno/hg19_gencode_from_ucsc_nothick_nocds_test.bed";
		String in_bed_path = "/Users/mengli/Documents/projects/iso/anno/hg19_gencode_from_ucsc_nothick_nocds.bed";

		String bam_path = "/Volumes/mengli/eric_nu/bam_3rd/SRR8932661_1.bam";
		
		String out_path = "/Users/mengli/Downloads/aaaa.txt";

		int intron_flank_threshold=20;
		
		iso(in_bed_path,out_path,bam_path,intron_flank_threshold);

		
		
		/*
		String in_bed_path = "/Users/mengli/Documents/projects/iso/zebrafish/GRCz11_ensembl_nothick.bed";

		String bam_path = "/Users/mengli/Downloads/";
		
		String out_path = "/Users/mengli/Documents/projects/iso/result/iso_test.tsv";

		int intron_flank_threshold=20;
		
		boolean unique_int=true;
		
		boolean consider_large_int=true;
		
		e_map=Reada5a3.read5a3("/Users/mengli/Documents/projects/iso/zebrafish/a3a5_bias_all_zebra.tsv");
		 */
		
	}
	
	public static void main(String[] args) throws Exception {	
		//test_iso();
		
        Options options = new Options();

        Option input = new Option("i", "bed", true, "annotation bed file path");
        input.setRequired(true);
        options.addOption(input);

        Option ibam = new Option("ibam", "inputbam", true, "input sorted and indexed bam file");
        ibam.setRequired(true);
        options.addOption(ibam);
        
        Option output = new Option("o", "output", true, "output file");
        output.setRequired(true);
        options.addOption(output);
        
        Option thredshold = new Option("t", "threshold", true, "the minimum number of base pair mapped into introns");
        thredshold.setRequired(false);
        options.addOption(thredshold);
    
       
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = null;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("intron splicing order pair calculation", options);

            System.exit(1);
        }

        String inputFilePath = cmd.getOptionValue("bed");
        
        String outputFilePath = cmd.getOptionValue("output");
        
        String input_bam_file = cmd.getOptionValue("inputbam");
        

        File input_file = new File(inputFilePath);
        if(!input_file.exists()) {
        	System.err.print("input bed annotation not exists");
            System.exit(1);
        }
        
       
        File output_file = new File(outputFilePath).getAbsoluteFile();
        if((!output_file.exists()) && (!output_file.getParentFile().exists())) {
        	System.err.print("output path not exitst");
            System.exit(1);
        }
        
        File bam_file = new File(input_bam_file);
        if(!bam_file.exists()) {
        	System.err.print("bam file not exitst");
            System.exit(1);
        }
        
        File bam_index_file = new File(input_bam_file+".bai");
        if(!bam_index_file.exists()) {
        	System.err.print("bam index file not exitst");
            System.exit(1);
        }
        
        
        int thredshold_v = 90;

        if( cmd.hasOption("threshold")  ) {
            String thredshold_str = cmd.getOptionValue("threshold");

        	if(Reada5a3.isInteger(thredshold_str) ){
        	
        		thredshold_v=Integer.parseInt( thredshold_str );
            
	        }else {
	        	System.err.print("threshold must be integer");
	
	            System.exit(1);
	
	        }
        }

        
        
		iso(inputFilePath,outputFilePath,input_bam_file,thredshold_v);
		
	}
	
	public static void iso(String in_bed_path,String out_path,String bam_path, int t_threshold  ) throws Exception {
		
		String bam_file_name=new File(bam_path).getName();
		
		System.out.println("bam file name: "+ bam_file_name);
		
		System.out.println("annotation bed file name: "+ new File(in_bed_path).getAbsolutePath() );
		
		System.out.println("output path name: "+ new File(out_path).getAbsolutePath());

		System.out.println("read map into intron threshold: "+ t_threshold+"bp");

		
		int total_intron_splicing_order_pair=0;
		
		HashMap<String, Integer> e_map;
		
		//String in_bed_path=t_input_path;
		//String bam_path=args_l[1];
		//String out_path = args_l[2];
		
		//if(args_l.length>3) {
			//intron_flank_threshold=Integer.parseInt(args_l[3]);
		//intron_flank_threshold=t_threshold;
		//}
		
		
		String out_path_read_names=bam_path+"_read_name.txt";
		
		//the minimum length of reads in anchor region
	
		
		//only consider unique intron
		boolean unique_int=false;
		
		//if consider all introns in transcript or just neighbor introns.
		boolean consider_large_int=true;
		
		//one isoform's exon maybe in another isoform's intron region
		boolean consider_exon_in_intron=true;

		//correct a5a3 bias, not used
		e_map=Reada5a3.read5a3_no5a3();
		
		
		//System.out.println("consider whole transcript: "+consider_large_int);
		//System.out.println("unique intron: "+unique_int);
				
		File bed_file = new File(in_bed_path);
				
		AbstractFeatureReader<BEDFeature, LineIterator> br = AbstractFeatureReader.getFeatureReader(bed_file.getPath(),
				new BEDCodec(), false);
		
		
		File bamFile = new File(bam_path);

		SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
        
		HashSet<String> contig_names=preprocess_name(reader);

		
		BufferedWriter output_one = new BufferedWriter(new FileWriter(out_path));

		BufferedWriter output_one_read_name = new BufferedWriter(new FileWriter(out_path_read_names) );

		
		HashMap<String, Integer> interest_iso = new HashMap<String, Integer>();

		HashMap<String, Integer> read_support_jc_pair = new HashMap<String, Integer>();

		
		Iterable<BEDFeature> iter = br.iterator();

		HashSet<String> two_int_coord = new HashSet<String>();
		
		HashSet<String> all_introns=Reada5a3.get_all_intron(bed_file);
		
        int bed_line_numb=0;
		for (BEDFeature feat : iter) {
			bed_line_numb++;
		}
		
        System.out.print(bed_line_numb);
        
		
        //String anim= "|/-\\";
        float x=0;
        
		iter = br.iterator();
		for (BEDFeature feat : iter) {
		
			String chr =feat.getChr();
			
			String data = "\r"+"BAM file: "+bam_file_name+". Percent of transcripts finished: "+ String.format("%.5f", x/bed_line_numb*100) +"%";
			
	        //String data = "\r" + anim.charAt( Math.round(x/bed_line_numb) % anim.length()) + " " + Math.round(x/bed_line_numb);
	        System.out.write(data.getBytes());
	        
			x++;
			
			// System.out.println(st);
			int count_of_exon = feat.getExons().size();
			if (count_of_exon < 3) {
				continue;
			}
			
			if(!feat.getChr().startsWith("chr") ) {
				continue;
			}
			
			if(!contig_names.contains(chr) ) {
				System.out.println("\nChromosome contig not found in BAM file: "+ chr);
				continue;
			}
			
			
			//if(feat.getChr().contains("gl") || feat.getChr().contains("Un_") ) {
			//	continue;
			//}
			
			String strand=feat.getStrand().toString();
			String trans_name=feat.getName();
			
			
			int trans_start=feat.getStart();
			
			int trans_end=feat.getEnd();
			

			/*
			 * record all intron coordiantes, if juntioin containing read doesn't compatible with 
			 * the intron, then they are not in this transcript. 
			 */
			HashSet<String> all_intron_trans=new HashSet<String>();
			
			for (int i = 1; i < (count_of_exon ); i++) {

				Exon e=feat.getExons().get(i);


				int start = e.getCdStart();

				int left_int_coor=feat.getExons().get(i-1).getCdEnd()+1;

				String left_int = chr+":"+left_int_coor+"-"+(start-1);
				
				all_intron_trans.add(left_int);
				
			}
			
			
			for (int i = 1; i < (count_of_exon ); i++) {

				Exon e=feat.getExons().get(i);
				
				//String chr =feat.getChr();
				
				int start = e.getCdStart();				
				
				int left_int_coor=feat.getExons().get(i-1).getCdEnd()+1;

				String left_int = chr+":"+left_int_coor+"-"+(start-1);
				
				
				HashSet<String> cover_left_names = new HashSet<String>();
				
				HashSet<String> left_jc_names = new HashSet<String>();
				
				
				HashSet<String> cover_left_names_intron = new HashSet<String>();

				HashSet<String> left_jc_names_intron = new HashSet<String>();

				
				if(unique_int && two_int_coord.contains(chr+":"+left_int) ) {
					continue;
				}
				
				two_int_coord.add(chr+":"+left_int);
				
				int right_exon_start=start;
				
				int left_exon_end=left_int_coor-1;
				
				boolean exon_ri_a3_a5_likely_bias=false;
				
				int anchor_region_len=100;

				int intron_flank_threshold= t_threshold;

				//if anchor_region_len longer than the intron
				if(anchor_region_len> (right_exon_start-left_exon_end-1)) {
					anchor_region_len=right_exon_start-left_exon_end-1;
				}
				
				if(intron_flank_threshold> (right_exon_start-left_exon_end-1)) {
					intron_flank_threshold=right_exon_start-left_exon_end-1-1;
				}
				
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
				
				
				CloseableIterator<SAMRecord> r = null;
				//System.out.println(chr);System.out.println(left_int_coor);System.out.println(start);

				try {
					r=reader.query(chr, left_int_coor-1, start+1,false);
				}catch(Exception e1) {
					System.out.println(e1.getMessage());
					e1.printStackTrace();
					System.out.println(chr);System.out.println(left_int_coor);System.out.println(start);

					continue;
				}
				
				
				while (r.hasNext()) {
					
					SAMRecord recode = r.next();
					
					//bwa, mapping quality >0 is unique map
					if( (recode.getMappingQuality()<1)  ) {
						continue;
					}
					
					if(recode.getReadLength()<30 ) {
						continue;
					}
					
					//star mapping, 255 is unique mapping
					//STAR has NH tag, >1 not unique, =1 unique mapping
					if(recode.hasAttribute("NH") ) {
						if(recode.getIntegerAttribute("NH")>1) {
							continue;
						}
						
					}
					
					String name = recode.getReadName();				
					
					boolean cover_left=false;
					//boolean cover_right=false;
					
					boolean has_N = recode.getCigar().containsOperator(CigarOperator.N);
					
					
					if(has_N) {
						
						int[] att =null;
						
						if(recode.hasAttribute("jI")) {
							att=recode.getSignedIntArrayAttribute("jI");
						}else {
							att=Reada5a3.get_splice_junction_pos(recode.getCigar(),recode.getStart() );
						}
						
						if(att.length<2) {
							continue;
						}
						
						
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
				        
				        //the minimum of intron length and intron flank threshold
				        //int intron_flank_threshold_adjust=Math.min(intron_flank_threshold,  start-left_int_coor );
				        
				        
				        // if the junctions in this read contain the intron, add this read to the left_jc_names
						for ( int j = 0; j < (Math.floor(att.length / 2) ); j++) {
							
							if((att[ j * 2] == left_int_coor) &&  (att[ j * 2+1] == start-1 ) ) {
				        		no_left_sj=false;
				        		
								if( (recode.getStart() < left_int_coor- intron_flank_threshold) &&
										(recode.getEnd() > start-1+intron_flank_threshold)	) {
									
					        		left_jc_names.add(name);
					        		
					        		break;
								}
							}
						}
				        
						if(no_left_sj==false) {
							cover_left=false;
							continue;
						}
				        
				        
				        //no junction site in the intron
				        for(int tt=0;tt<att.length;tt++) {
				        	if( (att[tt]<=start-1) && (att[tt]>=left_exon_end+1 ) ) {
				        		no_left_sj=false;
				        	}
				        }
				        
						if(no_left_sj==false) {
							cover_left=false;
							continue;
						}
				        
				        //no overlap junc span the intron,
				        //if no junc span the intron, either, the junction are in right or left of the intron
				        // if junc in both the right and left intron, then the read must cover the whole intron
				        if(no_left_sj==true) {
							for ( int j = 0; j < (Math.floor(att.length / 2) ); j++) {
	
								if((att[ j * 2]<=start-1) &&  (att[ j * 2+1]>=left_int_coor ) ) {
					        		no_left_sj=false;
					        		continue;
								}
							
							}
				        }
				        
						if(no_left_sj==false) {
							cover_left=false;
							continue;
						}
						
						
						if(consider_exon_in_intron) {
					        /*
					         *  force the read overlap with intron at least intron_flank_threshold bp. 
					         *  read must around with intron-exon junction ( -80bp->0) to ensure not spliced. 
					         *  and the minimum bases in this region is intron_flank_threshold
					         *  
					         */
								
							cover_left = (no_left_sj) && (
									( (recode.getStart() <= (right_exon_start - intron_flank_threshold)) &&
									(recode.getEnd() >= (right_exon_start - anchor_region_len+intron_flank_threshold) )    ) ||
									(   (recode.getStart() <=(left_exon_end + anchor_region_len-intron_flank_threshold) )  &&
									(recode.getEnd() >= (left_exon_end + intron_flank_threshold)) )
									);
							
						}else {
							cover_left = (no_left_sj) && (
									( (recode.getStart() <= (right_exon_start - intron_flank_threshold))    ) &&
									(   
									(recode.getEnd() >= (left_exon_end + intron_flank_threshold)) )
									);
						}
										
				 
					}else {
						if(consider_exon_in_intron) {

						cover_left = ( (recode.getStart() <= (right_exon_start - intron_flank_threshold)) &&
								(recode.getEnd() >= (right_exon_start - anchor_region_len+intron_flank_threshold) )    ) ||
								(   (recode.getStart() <= (left_exon_end + anchor_region_len-intron_flank_threshold) )  &&
								(recode.getEnd() >= (left_exon_end + intron_flank_threshold)) );
						}else {
							
							cover_left = ( (recode.getStart() <= (right_exon_start - intron_flank_threshold))  ) &&
									(  
									(recode.getEnd() >= (left_exon_end + intron_flank_threshold)) );
						}
				
						
					}
					
					if (cover_left) {
						cover_left_names.add(name);
					}
					
				}
				
				r.close();
				
				
				CloseableIterator<SAMRecord> r_sec = null;
				
				if(cover_left_names.size()==0) {
					//continue;
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
					
					
					if(recode.getReadLength()<30 ) {
						continue;
					}
					
					//STAR mapping, 255 is unique mapping
					//STAR has NH tag, >1 not unique, =1 unique mapping
					if(recode.hasAttribute("NH") ) {
						if(recode.getIntegerAttribute("NH")>1) {
							continue;
						}
					}
					
					String paired_name=recode.getReadName();
					
					if (cover_left_names.contains(paired_name) ) {
						//output_one_read_name.write(paired_name+"\n");
						
					}else if(!left_jc_names.contains(paired_name)) {
						continue;
					}
					
					
					int[] att =null;
					
					if(recode.hasAttribute("jI")) {
						att=recode.getSignedIntArrayAttribute("jI");
					}else {
						att=Reada5a3.get_splice_junction_pos(recode.getCigar(),recode.getStart() );
					}
					
					boolean has_jc = (att.length > 0);
					
					if (!has_jc) {
						
						continue;
					}
					
					//if no junctions break
					if (att.length < 2) {

						continue;
					}
					
					boolean is_all_junctions_in_intron=Reada5a3.all_junctions_in_intron(all_introns,all_intron_trans, chr, att);
			        
					if(!is_all_junctions_in_intron) {
						continue;
					}
								
					//String low_int = "";

					/*
					if (cover_left_names.contains(paired_name) ) {
						//low_int = left_int;
					} else {
						continue;
					}
					*/
					
					//output candidate read names for picard extract target reads
					//if the read don't related to this intron
					
					for ( int j = 0; j < (Math.floor(att.length / 2) ); j++) {
						
						String right_jc = recode.getContig() + ":" + (att[ j * 2] ) + "-"+ (att[ j * 2+1] );
						
						String one_region = trans_name+"\t"+left_int + "\t" + right_jc+"\t"+strand+"\t"+exon_ri_a3_a5_likely_bias;
						
						if(left_int == right_jc) {
							continue;
						}
						
						
						//if the read cover the intron
						if (cover_left_names.contains(paired_name) ) {
							
							//two read in same fragment may overlapped, the below sentence avoid this bias
							//the read may overlapped some introns but not some others.
							//For every read or read pair, only 
							if (cover_left_names_intron.contains(paired_name+"-"+right_jc) ) {
								continue;
							}else {
								cover_left_names_intron.add(paired_name+"-"+right_jc);
							}
							
							if (!interest_iso.containsKey(one_region)) {
								interest_iso.put(one_region, 1);
								continue;
							}
							
							int one_rcount = interest_iso.get(one_region);
							interest_iso.put(one_region, one_rcount + 1);
						}
						
						
						//if the read contain the junction which is just the intron
						if(left_jc_names.contains(paired_name) ) {
							
							//two read in same fragment may overlapped, the below sentence avoid this bias
							//the read may overlapped some introns but not some others.
							if (left_jc_names_intron.contains(paired_name+"-"+right_jc) ) {
								continue;
							}else {
								left_jc_names_intron.add(paired_name+"-"+right_jc);
							}
							
							//output_one_read_name.write(one_region+"\t"+paired_name+"\n");
							
							if (!read_support_jc_pair.containsKey(one_region)) {
								read_support_jc_pair.put(one_region, 1);
								continue;
							}
							
							int one_rcount = read_support_jc_pair.get(one_region);
							read_support_jc_pair.put(one_region, one_rcount + 1);
							
						}
							
					}
					

					

					
				}
				r_sec.close();
				
			}
			
			
			for (Entry<String, Integer> m : interest_iso.entrySet()) {
				
				int read_number_support_jc=0;
				if(read_support_jc_pair.containsKey(m.getKey())) {
					read_number_support_jc=read_support_jc_pair.get(m.getKey());
				}
				
				output_one.write(m.getKey() + "\t" + m.getValue()+"\t"+read_number_support_jc);
				
				output_one.newLine();
				
				total_intron_splicing_order_pair++;
				
			}
			
			interest_iso.clear();
			
			read_support_jc_pair.clear();
		}

		
		for (Entry<String, Integer> m : interest_iso.entrySet()) {
			
			int read_number_support_jc=0;
			if(read_support_jc_pair.containsKey(m.getKey())) {
				read_number_support_jc=read_support_jc_pair.get(m.getKey());
			}
			
			output_one.write(m.getKey() + "\t" + m.getValue()+"\t"+read_number_support_jc) ;
			
			output_one.newLine();
		}
		
		String data = "\r"+"BAM file: "+bam_file_name+". Percent of transcripts finished: "+ String.format("%.5f", 1.0*100.0) +"%";
		
		System.out.write(data.getBytes());
        
		System.out.println("\nBam file finished: "+bam_file_name+
				". Total intron splicing order pairs: "+ total_intron_splicing_order_pair );
		
		read_support_jc_pair.clear();

		interest_iso.clear();
		
		output_one.close();
		
		output_one_read_name.close();


	}
	
}
