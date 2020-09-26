package hmm;

import shared.Tools;

public class HMMSearchLine {

	public HMMSearchLine(byte[] line_){
		line=line_;
		
		//FORCE_JAVA_PARSE_DOUBLE=true //put this somewhere else

		int a=0, b=0;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 0: "+new String(line);
		name=new String(line, a, b-a);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 1: "+new String(line);
		field1=new String(line, a, b-a);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 2: "+new String(line);
		field2=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 3: "+new String(line);
		hmmName=new String(line, a, b-a);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 4: "+new String(line);
		accession=new String(line, a, b-a);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 5: "+new String(line);
		field5=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 6: "+new String(line);
		field6=Tools.parseDouble(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 7: "+new String(line);
		field7=Tools.parseFloat(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 8: "+new String(line);
		field8=Tools.parseFloat(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 9: "+new String(line);
		field9=Tools.parseFloat(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 10: "+new String(line);
		field10=Tools.parseFloat(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 11: "+new String(line);
		field11=Tools.parseDouble(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 12: "+new String(line);
		field12=Tools.parseDouble(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 13: "+new String(line);
		field13=Tools.parseFloat(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 14: "+new String(line);
		field14=Tools.parseFloat(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 15: "+new String(line);
		field15=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 16: "+new String(line);
		field16=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 17: "+new String(line);
		field17=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 18: "+new String(line);
		field18=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 19: "+new String(line);
		field19=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 20: "+new String(line);
		field20=Tools.parseInt(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 21: "+new String(line);
		field21=Tools.parseFloat(line, a, b);
		b++;
		a=b;
		
		while(b<line.length && line[b]!='\t'){b++;}
		assert(b>a) : "Missing field 22: "+new String(line);
		field22=new String(line, a, b-a);
		b++;
		a=b;
		
	}
	
	public final byte[] line;

	//protein_1            -            257 ATP-synt_A           PF00119.18   211   
	//1.9e-49  159.6  27.5   
	//1   1   7.3e-51   2.5e-49  159.2  27.5
	//3   210    41   250    38   251 0.87 -
	
	//field0 - protein_1
	String name;

	//field1 - -
	String field1;

	//field2 - 257
	int field2; //length?

	//field3 - ATP-synt_A
	String hmmName;

	//field4 - PF00119.18
	String accession;
	
	//field5 - 211 
	int field5;

	//field6 - 1.9e-49
	double field6;

	//field7 - 159.6
	float field7;

	//field8 - 27.5
	float field8;

	//field9 - 1
	float field9;

	//field10 - 1
	float field10;

	//field11 - 7.3e-51
	double field11;

	//field12 - 2.5e-49
	double field12;

	//field13 - 159.2
	float field13;

	//field14 - 27.5
	float field14;

	//field15 - 3
	int field15;

	//field16 - 210
	int field16;

	//field17 - 41
	int field17;
	
	//field18 - 250
	int field18;

	//field19 - 38
	int field19;

	//field20 - 251
	int field20;

	//field21 - 0.87
	float field21;

	//field22 - -
	String field22;
	
	
}
