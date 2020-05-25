.intel_syntax noprefix

// Registers that are used for parameter passing:
#define reg_p1  rdi
#define reg_p2  rsi
#define reg_p3  rdx
#define reg_p4  rcx

#define PARAM_BARR_DIV 30

#define LAZYR
#define BARRR

.text
//**********************************************************************************
//  Forward DGT, case N = 1024
//  Operation: c [reg_p1] <- NTT(a) [reg_p2], 
//             [reg_p3] points to table with the DGT constants. 
//  Output c and input a are in extended form. 
//********************************************************************************** 
.global poly_dgt_asm
poly_dgt_asm:
  xor          rax, rax 
  xor          r8, r8
  vmovdqa      ymm15, YMMWORD PTR [PARAM_Qx4+rip]
  vmovdqa      ymm14, YMMWORD PTR [PARAM_QINVx4+rip]
  vmovdqa      ymm13, YMMWORD PTR [PARAM_BARRx4+rip]
  mov          r10, 32                                
loop_dgt123:

// Round 1  
  vpmovsxdq    ymm0, XMMWORD PTR [reg_p2+rax]         // x[i] 
  vpmovsxdq    ymm1, XMMWORD PTR [reg_p2+rax+4*128]  
  vpmovsxdq    ymm2, XMMWORD PTR [reg_p2+rax+4*256]
  vpmovsxdq    ymm3, XMMWORD PTR [reg_p2+rax+4*384] 
  vpmovsxdq    ymm4, XMMWORD PTR [reg_p2+rax+4*512]    // x[i+m]           
  vpmovsxdq    ymm5, XMMWORD PTR [reg_p2+rax+4*640]      
  vpmovsxdq    ymm6, XMMWORD PTR [reg_p2+rax+4*768] 
  vpmovsxdq    ymm7, XMMWORD PTR [reg_p2+rax+4*896]
  
/*  vmovdqa    ymm0, YMMWORD PTR [reg_p2+2*rax]         // x[i] 
  vmovdqa    ymm1, YMMWORD PTR [reg_p2+2*rax+8*128]  
  vmovdqa    ymm2, YMMWORD PTR [reg_p2+2*rax+8*256]
  vmovdqa    ymm3, YMMWORD PTR [reg_p2+2*rax+8*384] 
  vmovdqa    ymm4, YMMWORD PTR [reg_p2+2*rax+8*512]    // x[i+m]           
  vmovdqa    ymm5, YMMWORD PTR [reg_p2+2*rax+8*640]      
  vmovdqa    ymm6, YMMWORD PTR [reg_p2+2*rax+8*768] 
  vmovdqa    ymm7, YMMWORD PTR [reg_p2+2*rax+8*896]       */
  vpsubq       ymm8, ymm0, ymm4 // sub = x[i] - x[i+m]
  vpsubq       ymm9, ymm1, ymm5
  vpsubq       ymm10, ymm2, ymm6
  vpsubq       ymm11, ymm3, ymm7
  vpaddq       ymm0, ymm0, ymm4 // add = (x[i] + x[i+m])
  vpaddq       ymm1, ymm1, ymm5
  vpaddq       ymm2, ymm2, ymm6
  vpaddq       ymm3, ymm3, ymm7
  
#ifndef BARRR
  vpsubd       ymm0, ymm0, ymm15                      // If result >= q then subtract q    
  vpsubd       ymm1, ymm1, ymm15                     
  vpsubd       ymm2, ymm2, ymm15                      
  vpsubd       ymm3, ymm3, ymm15    
  vpsrad       ymm4, ymm0, 31                 
  vpsrad       ymm5, ymm1, 31                       
  vpsrad       ymm6, ymm2, 31                       
  vpsrad       ymm7, ymm3, 31   
  vpand        ymm4, ymm4, ymm15                                             
  vpand        ymm5, ymm5, ymm15                       
  vpand        ymm6, ymm6, ymm15                   
  vpand        ymm7, ymm7, ymm15
  vpaddd       ymm0, ymm0, ymm4
  vpaddd       ymm1, ymm1, ymm5                     
  vpaddd       ymm2, ymm2, ymm6                      
  vpaddd       ymm3, ymm3, ymm7
#else
  vpmuldq      ymm4, ymm0, ymm13 // barr_reduce(add)
  vpmuldq      ymm5, ymm1, ymm13 // barr_reduce(add)
  vpmuldq      ymm6, ymm2, ymm13 // barr_reduce(add)
  vpmuldq      ymm7, ymm3, ymm13 // barr_reduce(add)
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuludq     ymm4, ymm4, ymm15
  vpmuludq     ymm5, ymm5, ymm15
  vpmuludq     ymm6, ymm6, ymm15
  vpmuludq     ymm7, ymm7, ymm15
  vpsubq       ymm0, ymm0, ymm4
  vpsubq       ymm1, ymm1, ymm5
  vpsubq       ymm2, ymm2, ymm6
  vpsubq       ymm3, ymm3, ymm7
#endif

  vpmovsxdq    ymm4, XMMWORD PTR [reg_p3+r8]
  vpmovsxdq    ymm5, XMMWORD PTR [reg_p3+r8+64*4]
  vpmovsxdq    ymm6, XMMWORD PTR [reg_p3+r8+128*4]
  vpmovsxdq    ymm7, XMMWORD PTR [reg_p3+r8+192*4]
  vpermq       ymm4, ymm4, 0x50
  vpermq       ymm5, ymm5, 0x50
  vpermq       ymm6, ymm6, 0x50
  vpermq       ymm7, ymm7, 0x50
  vpmuldq      ymm4, ymm8, ymm4
  vpmuldq      ymm5, ymm9, ymm5
  vpmuldq      ymm6, ymm10, ymm6
  vpmuldq      ymm7, ymm11, ymm7
  vpmuldq      ymm8, ymm4, ymm14 // reduce(a * sub)
  vpmuldq      ymm9, ymm5, ymm14 // reduce(a * sub)
  vpmuldq      ymm10, ymm6, ymm14 // reduce(a * sub)
  vpmuldq      ymm11, ymm7, ymm14 // reduce(a * sub)
  vpmuludq     ymm8, ymm8, ymm15
  vpmuludq     ymm9, ymm9, ymm15
  vpmuludq     ymm10, ymm10, ymm15
  vpmuludq     ymm11, ymm11, ymm15
  vpaddq       ymm4, ymm4, ymm8
  vpaddq       ymm5, ymm5, ymm9
  vpaddq       ymm6, ymm6, ymm10
  vpaddq       ymm7, ymm7, ymm11
  vpsrlq       ymm4, ymm4, 32
  vpsrlq       ymm5, ymm5, 32
  vpsrlq       ymm6, ymm6, 32
  vpsrlq       ymm7, ymm7, 32

//Round 2
  vpsubq       ymm8, ymm0, ymm2 // sub = x[i] - x[i+m]
  vpsubq       ymm9, ymm1, ymm3
  vpsubq       ymm10, ymm4, ymm6
  vpsubq       ymm11, ymm5, ymm7
  vpaddq       ymm0, ymm0, ymm2 // add = (x[i] + x[i+m])
  vpaddq       ymm1, ymm1, ymm3
  vpaddq       ymm2, ymm4, ymm6
  vpaddq       ymm3, ymm5, ymm7
  
#ifndef LAZYR
#ifndef BARRR
  vpsrad       ymm4, ymm0, 31                         // If result < 0 then add q                   
  vpsrad       ymm5, ymm1, 31                       
  vpsrad       ymm6, ymm2, 31                       
  vpsrad       ymm7, ymm3, 31   
  vpand        ymm4, ymm4, ymm15                                             
  vpand        ymm5, ymm5, ymm15                       
  vpand        ymm6, ymm6, ymm15                    
  vpand        ymm7, ymm7, ymm15
  vpaddd       ymm0, ymm4, ymm0                      
  vpaddd       ymm1, ymm5, ymm1                     
  vpaddd       ymm2, ymm6, ymm2
  vpaddd       ymm3, ymm7, ymm3
#else
  vpmuldq      ymm4, ymm0, ymm13 // barr_reduce(add)
  vpmuldq      ymm6, ymm1, ymm13
  vpmuldq      ymm5, ymm2, ymm13
  vpmuldq      ymm7, ymm3, ymm13
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuludq     ymm4, ymm4, ymm15
  vpmuludq     ymm5, ymm5, ymm15
  vpmuludq     ymm6, ymm6, ymm15
  vpmuludq     ymm7, ymm7, ymm15
  vpsubq       ymm0, ymm0, ymm4
  vpsubq       ymm1, ymm1, ymm6
  vpsubq       ymm2, ymm2, ymm5
  vpsubq       ymm3, ymm3, ymm7
#endif
#endif

  vpmovsxdq    ymm12, XMMWORD PTR [reg_p3+2*r8]
  vpunpcklqdq  ymm5, ymm12, ymm12
  vpmovsxdq    ymm12, XMMWORD PTR [reg_p3+2*r8+4*128]
  vpunpcklqdq  ymm12, ymm12, ymm12
  vpmuldq      ymm4, ymm8, ymm5
  vpmuldq      ymm5, ymm10, ymm5
  vpmuldq      ymm6, ymm9, ymm12
  vpmuldq      ymm7, ymm11, ymm12
  vpmuldq      ymm8, ymm4, ymm14 // reduce(a * sub)
  vpmuldq      ymm9, ymm5, ymm14
  vpmuldq      ymm10, ymm6, ymm14
  vpmuldq      ymm11, ymm7, ymm14
  vpmuludq     ymm8, ymm8, ymm15
  vpmuludq     ymm9, ymm9, ymm15
  vpmuludq     ymm10, ymm10, ymm15
  vpmuludq     ymm11, ymm11, ymm15
  vpaddq       ymm4, ymm4, ymm8
  vpaddq       ymm5, ymm5, ymm9
  vpaddq       ymm6, ymm6, ymm10
  vpaddq       ymm7, ymm7, ymm11
  vpsrlq       ymm4, ymm4, 32
  vpsrlq       ymm5, ymm5, 32
  vpsrlq       ymm6, ymm6, 32
  vpsrlq       ymm7, ymm7, 32
  
//Round 3
  vpsubq       ymm8, ymm0, ymm1 // sub = x[i] - x[i+m]
  vpsubq       ymm9, ymm4, ymm6
  vpsubq       ymm10, ymm2, ymm3
  vpsubq       ymm11, ymm5, ymm7
  vpaddq       ymm0, ymm0, ymm1 // add = (x[i] + x[i+m])
  vpaddq       ymm1, ymm4, ymm6
  vpaddq       ymm2, ymm2, ymm3
  vpaddq       ymm3, ymm5, ymm7

  vpmuldq      ymm4, ymm0, ymm13 // barr_reduce(add)
  vpmuldq      ymm5, ymm1, ymm13
  vpmuldq      ymm6, ymm2, ymm13
  vpmuldq      ymm7, ymm3, ymm13
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuludq     ymm4, ymm4, ymm15
  vpmuludq     ymm5, ymm5, ymm15
  vpmuludq     ymm6, ymm6, ymm15
  vpmuludq     ymm7, ymm7, ymm15
  vpsubq       ymm0, ymm0, ymm4
  vpsubq       ymm1, ymm1, ymm5
  vpsubq       ymm2, ymm2, ymm6
  vpsubq       ymm3, ymm3, ymm7

  vpmovsxdq    ymm12, XMMWORD PTR [reg_p3+4*r8]
  vpmovsxdq    ymm4, XMMWORD PTR [reg_p3+4*r8+16]
  vpunpcklqdq  ymm12, ymm12, ymm4
  vpermq       ymm12, ymm12, 0x50
  vpmuldq      ymm4, ymm8, ymm12
  vpmuldq      ymm5, ymm9, ymm12
  vpmuldq      ymm6, ymm10, ymm12
  vpmuldq      ymm7, ymm11, ymm12
  vpmuldq      ymm8, ymm4, ymm14 // reduce(a * sub)
  vpmuldq      ymm9, ymm5, ymm14
  vpmuldq      ymm10, ymm6, ymm14
  vpmuldq      ymm11, ymm7, ymm14
  vpmuludq     ymm8, ymm8, ymm15
  vpmuludq     ymm9, ymm9, ymm15
  vpmuludq     ymm10, ymm10, ymm15
  vpmuludq     ymm11, ymm11, ymm15
  vpaddq       ymm4, ymm4, ymm8
  vpaddq       ymm5, ymm5, ymm9
  vpaddq       ymm6, ymm6, ymm10
  vpaddq       ymm7, ymm7, ymm11
  vpsrlq       ymm4, ymm4, 32
  vpsrlq       ymm5, ymm5, 32
  vpsrlq       ymm6, ymm6, 32
  vpsrlq       ymm7, ymm7, 32
  
  vmovdqa      [reg_p1+2*rax],       ymm0
  vmovdqa      [reg_p1+2*rax+8*128], ymm4
  vmovdqa      [reg_p1+2*rax+8*256], ymm1
  vmovdqa      [reg_p1+2*rax+8*384], ymm5
  vmovdqa      [reg_p1+2*rax+8*512], ymm2
  vmovdqa      [reg_p1+2*rax+8*640], ymm6
  vmovdqa      [reg_p1+2*rax+8*768], ymm3
  vmovdqa      [reg_p1+2*rax+8*896], ymm7

  add          rax, 16
  add          r8, 8
  dec          r10
  jnz          loop_dgt123

  xor          rax, rax
  xor          r11, r11
  xor          r8, r8
round4_loop:
  mov          r10, 2
  vpmovsxdq    ymm12, XMMWORD PTR [reg_p3+4*r8]
  vpmovsxdq    ymm4, XMMWORD PTR [reg_p3+4*r8+4*8]
  vpunpcklqdq  ymm12, ymm12, ymm4
  vpermq       ymm12, ymm12, 0x50
loop_dgt4:
  vmovdqa      ymm0, [reg_p1+rax]
  vmovdqa      ymm1, [reg_p1+rax+8*64]
  vmovdqa      ymm2, [reg_p1+rax+8*128]
  vmovdqa      ymm3, [reg_p1+rax+8*192]
  vmovdqa      ymm4, [reg_p1+rax+8*256]
  vmovdqa      ymm5, [reg_p1+rax+8*320]
  vmovdqa      ymm6, [reg_p1+rax+8*384]
  vmovdqa      ymm7, [reg_p1+rax+8*448]
  vpsubq       ymm8, ymm0, ymm1 // sub = x[i] - x[i+m]
  vpsubq       ymm9, ymm2, ymm3
  vpsubq       ymm10, ymm4, ymm5
  vpsubq       ymm11, ymm6, ymm7
  vpaddq       ymm0, ymm0, ymm1 // add = (x[i] + x[i+m])
  vpaddq       ymm1, ymm2, ymm3
  vpaddq       ymm2, ymm4, ymm5
  vpaddq       ymm3, ymm6, ymm7
  
#ifndef LAZYR
#ifndef BARRR
  vpsrad       ymm4, ymm0, 31                         // If result < 0 then add q                   
  vpsrad       ymm5, ymm1, 31                       
  vpsrad       ymm6, ymm2, 31                       
  vpsrad       ymm7, ymm3, 31   
  vpand        ymm4, ymm4, ymm15                                             
  vpand        ymm5, ymm5, ymm15                       
  vpand        ymm6, ymm6, ymm15                    
  vpand        ymm7, ymm7, ymm15
  vpaddd       ymm0, ymm4, ymm0                      
  vpaddd       ymm1, ymm5, ymm1                     
  vpaddd       ymm2, ymm6, ymm2
  vpaddd       ymm3, ymm7, ymm3
#else
  vpmuldq      ymm4, ymm0, ymm13 // barr_reduce(add)
  vpmuldq      ymm5, ymm1, ymm13 // barr_reduce(add)
  vpmuldq      ymm6, ymm2, ymm13 // barr_reduce(add)
  vpmuldq      ymm7, ymm3, ymm13 // barr_reduce(add)
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuludq     ymm4, ymm4, ymm15
  vpmuludq     ymm5, ymm5, ymm15
  vpmuludq     ymm6, ymm6, ymm15
  vpmuludq     ymm7, ymm7, ymm15
  vpsubq       ymm0, ymm0, ymm4
  vpsubq       ymm1, ymm1, ymm5
  vpsubq       ymm2, ymm2, ymm6
  vpsubq       ymm3, ymm3, ymm7
#endif
#endif

  vpmuldq      ymm4, ymm8, ymm12
  vpmuldq      ymm5, ymm9, ymm12
  vpmuldq      ymm6, ymm10, ymm12
  vpmuldq      ymm7, ymm11, ymm12
  vpmuldq      ymm8, ymm4, ymm14 // reduce(a * sub)
  vpmuldq      ymm9, ymm5, ymm14
  vpmuldq      ymm10, ymm6, ymm14
  vpmuldq      ymm11, ymm7, ymm14
  vpmuludq     ymm8, ymm8, ymm15
  vpmuludq     ymm9, ymm9, ymm15
  vpmuludq     ymm10, ymm10, ymm15
  vpmuludq     ymm11, ymm11, ymm15
  vpaddq       ymm4, ymm4, ymm8
  vpaddq       ymm5, ymm5, ymm9
  vpaddq       ymm6, ymm6, ymm10
  vpaddq       ymm7, ymm7, ymm11
  vpsrlq       ymm4, ymm4, 32
  vpsrlq       ymm5, ymm5, 32
  vpsrlq       ymm6, ymm6, 32
  vpsrlq       ymm7, ymm7, 32
  
  vmovdqa      [reg_p1+rax],       ymm0
  vmovdqa      [reg_p1+rax+8*64],  ymm4
  vmovdqa      [reg_p1+rax+8*128], ymm1
  vmovdqa      [reg_p1+rax+8*192], ymm5
  vmovdqa      [reg_p1+rax+8*256], ymm2
  vmovdqa      [reg_p1+rax+8*320], ymm6
  vmovdqa      [reg_p1+rax+8*384], ymm3
  vmovdqa      [reg_p1+rax+8*448], ymm7

  add          rax, 8*512
  dec          r10
  jnz          loop_dgt4
  add          r8, 16
  add          r11, 4*8
  mov          rax, r11
  cmp          r11, 16*(4*8)
  jne          round4_loop
  
  xor          rax, rax
  xor          r11, r11
  xor          r8, r8
round5_loop:
  mov          r10, 4
  vpmovsxdq    ymm12, XMMWORD PTR [reg_p3+4*r8]
  vpmovsxdq    ymm4, XMMWORD PTR [reg_p3+4*r8+4*16]
  vpunpcklqdq  ymm12, ymm12, ymm4
  vpermq       ymm12, ymm12, 0x50
loop_dgt5:
  vmovdqa      ymm0, [reg_p1+rax]
  vmovdqa      ymm1, [reg_p1+rax+8*32]
  vmovdqa      ymm2, [reg_p1+rax+8*64]
  vmovdqa      ymm3, [reg_p1+rax+8*96]
  vmovdqa      ymm4, [reg_p1+rax+8*128]
  vmovdqa      ymm5, [reg_p1+rax+8*160]
  vmovdqa      ymm6, [reg_p1+rax+8*192]
  vmovdqa      ymm7, [reg_p1+rax+8*224]
  vpsubq       ymm8, ymm0, ymm1 // sub = x[i] - x[i+m]
  vpsubq       ymm9, ymm2, ymm3
  vpsubq       ymm10, ymm4, ymm5
  vpsubq       ymm11, ymm6, ymm7
  vpaddq       ymm0, ymm0, ymm1 // add = (x[i] + x[i+m])
  vpaddq       ymm1, ymm2, ymm3
  vpaddq       ymm2, ymm4, ymm5
  vpaddq       ymm3, ymm6, ymm7
  
  vpmuldq      ymm4, ymm0, ymm13 // barr_reduce(add)
  vpmuldq      ymm5, ymm1, ymm13 // barr_reduce(add)
  vpmuldq      ymm6, ymm2, ymm13 // barr_reduce(add)
  vpmuldq      ymm7, ymm3, ymm13 // barr_reduce(add)
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuludq     ymm4, ymm4, ymm15
  vpmuludq     ymm5, ymm5, ymm15
  vpmuludq     ymm6, ymm6, ymm15
  vpmuludq     ymm7, ymm7, ymm15
  vpsubq       ymm0, ymm0, ymm4
  vpsubq       ymm1, ymm1, ymm5
  vpsubq       ymm2, ymm2, ymm6
  vpsubq       ymm3, ymm3, ymm7
  
  vpmuldq      ymm4, ymm8, ymm12
  vpmuldq      ymm5, ymm9, ymm12
  vpmuldq      ymm6, ymm10, ymm12
  vpmuldq      ymm7, ymm11, ymm12
  vpmuldq      ymm8, ymm4, ymm14 // reduce(a * sub)
  vpmuldq      ymm9, ymm5, ymm14
  vpmuldq      ymm10, ymm6, ymm14
  vpmuldq      ymm11, ymm7, ymm14
  vpmuludq     ymm8, ymm8, ymm15
  vpmuludq     ymm9, ymm9, ymm15
  vpmuludq     ymm10, ymm10, ymm15
  vpmuludq     ymm11, ymm11, ymm15
  vpaddq       ymm4, ymm4, ymm8
  vpaddq       ymm5, ymm5, ymm9
  vpaddq       ymm6, ymm6, ymm10
  vpaddq       ymm7, ymm7, ymm11
  vpsrlq       ymm4, ymm4, 32
  vpsrlq       ymm5, ymm5, 32
  vpsrlq       ymm6, ymm6, 32
  vpsrlq       ymm7, ymm7, 32
  
  vmovdqa      [reg_p1+rax],       ymm0
  vmovdqa      [reg_p1+rax+8*32],  ymm4
  vmovdqa      [reg_p1+rax+8*64], ymm1
  vmovdqa      [reg_p1+rax+8*96], ymm5
  vmovdqa      [reg_p1+rax+8*128], ymm2
  vmovdqa      [reg_p1+rax+8*160], ymm6
  vmovdqa      [reg_p1+rax+8*192], ymm3
  vmovdqa      [reg_p1+rax+8*224], ymm7

  add          rax, 8*256
  dec          r10
  jnz          loop_dgt5
  add          r8, 32
  add          r11, 4*8
  mov          rax, r11
  cmp          r11, 8*(4*8)
  jne          round5_loop

  xor          rax, rax
  mov          r10, 32
loop_dgt6789:
// Round 6
  vmovdqa      ymm0, [reg_p1+2*rax]
  vmovdqa      ymm1, [reg_p1+2*rax+8*4]
  vmovdqa      ymm2, [reg_p1+2*rax+8*8]
  vmovdqa      ymm3, [reg_p1+2*rax+8*12]
  vmovdqa      ymm4, [reg_p1+2*rax+8*16]
  vmovdqa      ymm5, [reg_p1+2*rax+8*20]
  vmovdqa      ymm6, [reg_p1+2*rax+8*24]
  vmovdqa      ymm7, [reg_p1+2*rax+8*28]
  vpsubq       ymm8, ymm0, ymm4 // sub = x[i] - x[i+m]
  vpsubq       ymm9, ymm1, ymm5
  vpsubq       ymm10, ymm2, ymm6
  vpsubq       ymm11, ymm3, ymm7
  vpaddq       ymm0, ymm0, ymm4 // add = (x[i] + x[i+m])
  vpaddq       ymm1, ymm1, ymm5
  vpaddq       ymm2, ymm2, ymm6
  vpaddq       ymm3, ymm3, ymm7

#ifndef LAZYR
#ifndef BARRR
  vpsrad       ymm4, ymm0, 31                         // If result < 0 then add q                   
  vpsrad       ymm5, ymm1, 31                       
  vpsrad       ymm6, ymm2, 31                       
  vpsrad       ymm7, ymm3, 31   
  vpand        ymm4, ymm4, ymm15                                             
  vpand        ymm5, ymm5, ymm15                       
  vpand        ymm6, ymm6, ymm15                    
  vpand        ymm7, ymm7, ymm15
  vpaddd       ymm0, ymm4, ymm0                      
  vpaddd       ymm1, ymm5, ymm1                     
  vpaddd       ymm2, ymm6, ymm2
  vpaddd       ymm3, ymm7, ymm3
#else
  vpmuldq      ymm4, ymm0, ymm13 // barr_reduce(add)
  vpmuldq      ymm5, ymm1, ymm13 // barr_reduce(add)
  vpmuldq      ymm6, ymm2, ymm13 // barr_reduce(add)
  vpmuldq      ymm7, ymm3, ymm13 // barr_reduce(add)
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuludq     ymm4, ymm4, ymm15
  vpmuludq     ymm5, ymm5, ymm15
  vpmuludq     ymm6, ymm6, ymm15
  vpmuludq     ymm7, ymm7, ymm15
  vpsubq       ymm0, ymm0, ymm4
  vpsubq       ymm1, ymm1, ymm5
  vpsubq       ymm2, ymm2, ymm6
  vpsubq       ymm3, ymm3, ymm7
#endif
#endif
  
  vpmovsxdq    ymm12, XMMWORD PTR [reg_p3]
  vpmovsxdq    ymm4, XMMWORD PTR [reg_p3+128]
  vpunpcklqdq  ymm12, ymm12, ymm4
  vpermq       ymm12, ymm12, 0x50
  vpmuldq      ymm4, ymm8, ymm12
  vpmovsxdq    ymm12, XMMWORD PTR [reg_p3+256]
  vpmovsxdq    ymm5, XMMWORD PTR [reg_p3+384]
  vpunpcklqdq  ymm12, ymm12, ymm5
  vpermq       ymm12, ymm12, 0x50
  vpmuldq      ymm5, ymm9, ymm12
  vpmovsxdq    ymm12, XMMWORD PTR [reg_p3+512]
  vpmovsxdq    ymm6, XMMWORD PTR [reg_p3+640]
  vpunpcklqdq  ymm12, ymm12, ymm6
  vpermq       ymm12, ymm12, 0x50
  vpmuldq      ymm6, ymm10, ymm12
  vpmovsxdq    ymm12, XMMWORD PTR [reg_p3+768]
  vpmovsxdq    ymm7, XMMWORD PTR [reg_p3+896]
  vpunpcklqdq  ymm12, ymm12, ymm7
  vpermq       ymm12, ymm12, 0x50
  vpmuldq      ymm7, ymm11, ymm12
  vpmuldq      ymm8, ymm4, ymm14 // reduce(a * sub)
  vpmuldq      ymm9, ymm5, ymm14
  vpmuldq      ymm10, ymm6, ymm14
  vpmuldq      ymm11, ymm7, ymm14
  vpmuludq     ymm8, ymm8, ymm15
  vpmuludq     ymm9, ymm9, ymm15
  vpmuludq     ymm10, ymm10, ymm15
  vpmuludq     ymm11, ymm11, ymm15
  vpaddq       ymm4, ymm4, ymm8
  vpaddq       ymm5, ymm5, ymm9
  vpaddq       ymm6, ymm6, ymm10
  vpaddq       ymm7, ymm7, ymm11
  vpsrlq       ymm4, ymm4, 32
  vpsrlq       ymm5, ymm5, 32
  vpsrlq       ymm6, ymm6, 32
  vpsrlq       ymm7, ymm7, 32
  
// Round 7
  vpsubq       ymm8, ymm0, ymm2 // sub = x[i] - x[i+m]
  vpsubq       ymm9, ymm1, ymm3
  vpsubq       ymm10, ymm4, ymm6
  vpsubq       ymm11, ymm5, ymm7
  vpaddq       ymm0, ymm0, ymm2 // add = (x[i] + x[i+m])
  vpaddq       ymm1, ymm1, ymm3
  vpaddq       ymm2, ymm4, ymm6
  vpaddq       ymm3, ymm5, ymm7
  
  vpmuldq      ymm4, ymm0, ymm13 // barr_reduce(add)
  vpmuldq      ymm5, ymm1, ymm13 // barr_reduce(add)
  vpmuldq      ymm6, ymm2, ymm13 // barr_reduce(add)
  vpmuldq      ymm7, ymm3, ymm13 // barr_reduce(add)
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuludq     ymm4, ymm4, ymm15
  vpmuludq     ymm5, ymm5, ymm15
  vpmuludq     ymm6, ymm6, ymm15
  vpmuludq     ymm7, ymm7, ymm15
  vpsubq       ymm0, ymm0, ymm4
  vpsubq       ymm1, ymm1, ymm5
  vpsubq       ymm2, ymm2, ymm6
  vpsubq       ymm3, ymm3, ymm7
  
  vpmovsxdq    ymm12, XMMWORD PTR [reg_p3]
  vpmovsxdq    ymm4, XMMWORD PTR [reg_p3+256]
  vpunpcklqdq  ymm12, ymm12, ymm4
  vpermq       ymm12, ymm12, 0x50
  vpmuldq      ymm4, ymm8, ymm12
  vpmuldq      ymm6, ymm10, ymm12
  vpmovsxdq    ymm12, XMMWORD PTR [reg_p3+512]
  vpmovsxdq    ymm5, XMMWORD PTR [reg_p3+768]
  vpunpcklqdq  ymm12, ymm12, ymm5
  vpermq       ymm12, ymm12, 0x50
  vpmuldq      ymm5, ymm9, ymm12
  vpmuldq      ymm7, ymm11, ymm12
  vpmuldq      ymm8, ymm4, ymm14 // reduce(a * sub)
  vpmuldq      ymm9, ymm5, ymm14
  vpmuldq      ymm10, ymm6, ymm14
  vpmuldq      ymm11, ymm7, ymm14
  vpmuludq     ymm8, ymm8, ymm15
  vpmuludq     ymm9, ymm9, ymm15
  vpmuludq     ymm10, ymm10, ymm15
  vpmuludq     ymm11, ymm11, ymm15
  vpaddq       ymm4, ymm4, ymm8
  vpaddq       ymm5, ymm5, ymm9
  vpaddq       ymm6, ymm6, ymm10
  vpaddq       ymm7, ymm7, ymm11
  vpsrlq       ymm4, ymm4, 32
  vpsrlq       ymm5, ymm5, 32
  vpsrlq       ymm6, ymm6, 32
  vpsrlq       ymm7, ymm7, 32
  
// Round 8
  vpsubq       ymm8, ymm0, ymm1 // sub = x[i] - x[i+m]
  vpsubq       ymm9, ymm4, ymm5
  vpsubq       ymm10, ymm2, ymm3
  vpsubq       ymm11, ymm6, ymm7
  vpaddq       ymm0, ymm0, ymm1 // add = (x[i] + x[i+m])
  vpaddq       ymm1, ymm4, ymm5
  vpaddq       ymm2, ymm2, ymm3
  vpaddq       ymm3, ymm6, ymm7
  
#ifndef LAZYR
#ifndef BARRR
  vpsrad       ymm4, ymm0, 31                         // If result < 0 then add q                   
  vpsrad       ymm5, ymm1, 31                       
  vpsrad       ymm6, ymm2, 31                       
  vpsrad       ymm7, ymm3, 31   
  vpand        ymm4, ymm4, ymm15                                             
  vpand        ymm5, ymm5, ymm15                       
  vpand        ymm6, ymm6, ymm15                    
  vpand        ymm7, ymm7, ymm15
  vpaddd       ymm0, ymm4, ymm0                      
  vpaddd       ymm1, ymm5, ymm1                     
  vpaddd       ymm2, ymm6, ymm2
  vpaddd       ymm3, ymm7, ymm3
#else
  vpmuldq      ymm4, ymm0, ymm13 // barr_reduce(add)
  vpmuldq      ymm5, ymm1, ymm13 // barr_reduce(add)
  vpmuldq      ymm6, ymm2, ymm13 // barr_reduce(add)
  vpmuldq      ymm7, ymm3, ymm13 // barr_reduce(add)
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuludq     ymm4, ymm4, ymm15
  vpmuludq     ymm5, ymm5, ymm15
  vpmuludq     ymm6, ymm6, ymm15
  vpmuludq     ymm7, ymm7, ymm15
  vpsubq       ymm0, ymm0, ymm4
  vpsubq       ymm1, ymm1, ymm5
  vpsubq       ymm2, ymm2, ymm6
  vpsubq       ymm3, ymm3, ymm7
#endif
#endif
  
  vpmovsxdq    ymm12, XMMWORD PTR [reg_p3]
  vpmovsxdq    ymm4, XMMWORD PTR [reg_p3+512]
  vpunpcklqdq  ymm12, ymm12, ymm4
  vpermq       ymm12, ymm12, 0x50
  vpmuldq      ymm4, ymm8, ymm12
  vpmuldq      ymm6, ymm10, ymm12
  vpmuldq      ymm5, ymm9, ymm12
  vpmuldq      ymm7, ymm11, ymm12
  vpmuldq      ymm8, ymm4, ymm14 // reduce(a * sub)
  vpmuldq      ymm9, ymm5, ymm14
  vpmuldq      ymm10, ymm6, ymm14
  vpmuldq      ymm11, ymm7, ymm14
  vpmuludq     ymm8, ymm8, ymm15
  vpmuludq     ymm9, ymm9, ymm15
  vpmuludq     ymm10, ymm10, ymm15
  vpmuludq     ymm11, ymm11, ymm15
  vpaddq       ymm4, ymm4, ymm8
  vpaddq       ymm5, ymm5, ymm9
  vpaddq       ymm6, ymm6, ymm10
  vpaddq       ymm7, ymm7, ymm11
  vpsrlq       ymm4, ymm4, 32
  vpsrlq       ymm5, ymm5, 32
  vpsrlq       ymm6, ymm6, 32
  vpsrlq       ymm7, ymm7, 32
  
// Round 9
  vperm2i128   ymm12, ymm0, ymm4, 0x20              
  vperm2i128   ymm0, ymm0, ymm4, 0x31               
  vperm2i128   ymm4, ymm1, ymm5, 0x20
  vperm2i128   ymm1, ymm1, ymm5, 0x31
  vperm2i128   ymm5, ymm2, ymm6, 0x20
  vperm2i128   ymm2, ymm2, ymm6, 0x31
  vperm2i128   ymm6, ymm3, ymm7, 0x20
  vperm2i128   ymm3, ymm3, ymm7, 0x31
  vpsubq       ymm8, ymm12, ymm0 // sub = x[i] - x[i+m]
  vpsubq       ymm9, ymm4, ymm1
  vpsubq       ymm10, ymm5, ymm2
  vpsubq       ymm11, ymm6, ymm3
  vpaddq       ymm0, ymm12, ymm0 // add = (x[i] + x[i+m])
  vpaddq       ymm1, ymm4, ymm1
  vpaddq       ymm2, ymm5, ymm2
  vpaddq       ymm3, ymm6, ymm3
  
  vbroadcastss ymm12, DWORD PTR [reg_p3]
  
  vpmuldq      ymm4, ymm0, ymm13 // barr_reduce(add)
  vpmuldq      ymm5, ymm1, ymm13 // barr_reduce(add)
  vpmuldq      ymm6, ymm2, ymm13 // barr_reduce(add)
  vpmuldq      ymm7, ymm3, ymm13 // barr_reduce(add)
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuludq     ymm4, ymm4, ymm15
  vpmuludq     ymm5, ymm5, ymm15
  vpmuludq     ymm6, ymm6, ymm15
  vpmuludq     ymm7, ymm7, ymm15
  vpsubq       ymm0, ymm0, ymm4
  vpsubq       ymm1, ymm1, ymm5
  vpsubq       ymm2, ymm2, ymm6
  vpsubq       ymm3, ymm3, ymm7
  
  vpmuldq      ymm4, ymm8, ymm12
  vpmuldq      ymm5, ymm9, ymm12
  vpmuldq      ymm6, ymm10, ymm12
  vpmuldq      ymm7, ymm11, ymm12
  vpmuldq      ymm8, ymm4, ymm14 // reduce(a * sub)
  vpmuldq      ymm9, ymm5, ymm14
  vpmuldq      ymm10, ymm6, ymm14
  vpmuldq      ymm11, ymm7, ymm14
  vpmuludq     ymm8, ymm8, ymm15
  vpmuludq     ymm9, ymm9, ymm15
  vpmuludq     ymm10, ymm10, ymm15
  vpmuludq     ymm11, ymm11, ymm15
  vpaddq       ymm4, ymm4, ymm8
  vpaddq       ymm5, ymm5, ymm9
  vpaddq       ymm6, ymm6, ymm10
  vpaddq       ymm7, ymm7, ymm11
  vpsrlq       ymm4, ymm4, 32
  vpsrlq       ymm5, ymm5, 32
  vpsrlq       ymm6, ymm6, 32
  vpsrlq       ymm7, ymm7, 32
  
  vperm2i128   ymm12, ymm0, ymm4, 0x20              
  vperm2i128   ymm0, ymm0, ymm4, 0x31               
  vperm2i128   ymm4, ymm1, ymm5, 0x20
  vperm2i128   ymm1, ymm1, ymm5, 0x31
  vperm2i128   ymm5, ymm2, ymm6, 0x20
  vperm2i128   ymm2, ymm2, ymm6, 0x31
  vperm2i128   ymm6, ymm3, ymm7, 0x20
  vperm2i128   ymm3, ymm3, ymm7, 0x31
  
  vmovdqa      [reg_p1+2*rax],     ymm12
  vmovdqa      [reg_p1+2*rax+8*4], ymm0
  vmovdqa      [reg_p1+2*rax+8*8], ymm4
  vmovdqa      [reg_p1+2*rax+8*12], ymm1
  vmovdqa      [reg_p1+2*rax+8*16], ymm5
  vmovdqa      [reg_p1+2*rax+8*20], ymm2
  vmovdqa      [reg_p1+2*rax+8*24], ymm6
  vmovdqa      [reg_p1+2*rax+8*28], ymm3

  add          rax, 128
  dec          r10
  jnz          loop_dgt6789
  ret
  
//**********************************************************************************
//  Pointwise multiplication, case N = 1024
//  Operation: c [reg_p1] <- a [reg_p2] x b [reg_p3].
//  Output c is in extended form. 
//********************************************************************************** 
.global poly_pmul_asm
poly_pmul_asm:
  xor          rax, rax 
  vmovdqa      ymm15, YMMWORD PTR [PARAM_Qx4+rip]
  vmovdqa      ymm14, YMMWORD PTR [PARAM_QINVx4+rip]
  vmovdqa      ymm13, YMMWORD PTR [PERM_MASK+rip]
  mov          r10, 64                               
loop_pmul:
  vpmovsxdq      ymm0, XMMWORD PTR [reg_p2+rax]                     // Load a[j]
  vpmovsxdq      ymm1, XMMWORD PTR [reg_p2+rax+4*4]
  vpmovsxdq      ymm2, XMMWORD PTR [reg_p2+rax+4*512]
  vpmovsxdq      ymm3, XMMWORD PTR [reg_p2+rax+4*512+4*4]
  vpmovsxdq      ymm4, XMMWORD PTR [reg_p3+rax]                      
  vpmovsxdq      ymm5, XMMWORD PTR [reg_p3+rax+4*4]
  vpmovsxdq      ymm6, XMMWORD PTR [reg_p3+rax+4*512]
  vpmovsxdq      ymm7, XMMWORD PTR [reg_p3+rax+4*512+4*4]
  
  vpmuldq      ymm8, ymm0, ymm4
  vpmuldq      ymm9, ymm1, ymm5
  vpmuldq      ymm10, ymm2, ymm6
  vpmuldq      ymm11, ymm3, ymm7
  vpmuldq      ymm0, ymm0, ymm6
  vpmuldq      ymm1, ymm1, ymm7
  vpmuldq      ymm2, ymm2, ymm4
  vpmuldq      ymm3, ymm3, ymm5
  vpsubq       ymm8, ymm8, ymm10
  vpsubq       ymm9, ymm9, ymm11
  vpaddq       ymm10, ymm0, ymm2
  vpaddq       ymm11, ymm1, ymm3
  vpermq       ymm8, ymm8, 0b11011000
  vpermq       ymm10, ymm10, 0b11011000
  vpermq       ymm9, ymm9, 0b11011000
  vpermq       ymm11, ymm11, 0b11011000
  vpunpcklqdq     ymm12, ymm8, ymm10
  vpunpckhqdq     ymm8, ymm8, ymm10
  vpunpcklqdq     ymm10, ymm9, ymm11  
  vpunpckhqdq     ymm9, ymm9, ymm11  

  vpmuldq      ymm4, ymm12, ymm14 // reduce(a * sub)
  vpmuldq      ymm5, ymm10, ymm14
  vpmuldq      ymm6, ymm8, ymm14
  vpmuldq      ymm7, ymm9, ymm14
  vpmuludq     ymm4, ymm4, ymm15
  vpmuludq     ymm5, ymm5, ymm15
  vpmuludq     ymm6, ymm6, ymm15
  vpmuludq     ymm7, ymm7, ymm15
  vpaddq       ymm4, ymm4, ymm12
  vpaddq       ymm5, ymm5, ymm10
  vpaddq       ymm6, ymm6, ymm8
  vpaddq       ymm7, ymm7, ymm9
  vpermd       ymm4, ymm13, ymm4
  vpermd       ymm5, ymm13, ymm5
  vpermd       ymm6, ymm13, ymm6
  vpermd       ymm7, ymm13, ymm7
  vmovdqa      XMMWORD PTR [reg_p1+2*rax], xmm4
  vmovdqa      XMMWORD PTR [reg_p1+2*rax+4*4], xmm6
  vmovdqa      XMMWORD PTR [reg_p1+2*rax+4*8], xmm5
  vmovdqa      XMMWORD PTR [reg_p1+2*rax+4*12], xmm7

  add          rax, 4*8
  dec          r10
  jnz          loop_pmul
  ret

//**********************************************************************************
//  Pointwise multiplication, case N = 1024
//  Operation: c [reg_p1] <- a [reg_p2] x b [reg_p3].
//  Input b is in extended form. 
//********************************************************************************** 
.global poly_pmul_asm2
poly_pmul_asm2:
  xor          rax, rax
  vmovdqa      ymm15, YMMWORD PTR [PARAM_Qx4+rip]
  vmovdqa      ymm14, YMMWORD PTR [PARAM_QINVx4+rip]
  vmovdqa      ymm12, YMMWORD PTR [PERM_MASK+rip]
  mov          r10, 64                              
loop_pmul2:
  vpmovsxdq      ymm0, XMMWORD PTR [reg_p2+rax]                     // Load a[j]
  vpmovsxdq      ymm1, XMMWORD PTR [reg_p2+rax+4*4]
  vpmovsxdq      ymm2, XMMWORD PTR [reg_p2+rax+4*512]
  vpmovsxdq      ymm3, XMMWORD PTR [reg_p2+rax+4*512+4*4]
  vmovdqa      ymm4, YMMWORD PTR [reg_p3+4*rax]                      
  vmovdqa      ymm5, YMMWORD PTR [reg_p3+4*rax+4*8]
  vmovdqa      ymm6, YMMWORD PTR [reg_p3+4*rax+4*16]
  vmovdqa      ymm7, YMMWORD PTR [reg_p3+4*rax+4*24]
  vpunpcklqdq     ymm8, ymm4, ymm5
  vpunpckhqdq     ymm9, ymm4, ymm5
  vpunpcklqdq     ymm10, ymm6, ymm7
  vpunpckhqdq     ymm11, ymm6, ymm7
  vpermq       ymm8, ymm8, 0b11011000
  vpermq       ymm9, ymm9, 0b11011000
  vpermq       ymm10, ymm10, 0b11011000
  vpermq       ymm11, ymm11, 0b11011000
  
  vpmuldq      ymm4, ymm0, ymm8
  vpmuldq      ymm5, ymm1, ymm10
  vpmuldq      ymm6, ymm2, ymm9
  vpmuldq      ymm7, ymm3, ymm11
  vpmuldq      ymm0, ymm0, ymm9
  vpmuldq      ymm1, ymm1, ymm11
  vpmuldq      ymm2, ymm2, ymm8
  vpmuldq      ymm3, ymm3, ymm10
  vpsubq       ymm8, ymm4, ymm6
  vpsubq       ymm9, ymm5, ymm7
  vpaddq       ymm10, ymm0, ymm2
  vpaddq       ymm11, ymm1, ymm3
  vpermq       ymm8, ymm8, 0b11011000
  vpermq       ymm9, ymm9, 0b11011000
  vpermq       ymm10, ymm10, 0b11011000
  vpermq       ymm11, ymm11, 0b11011000
  vpunpcklqdq     ymm0, ymm8, ymm10
  vpunpckhqdq     ymm1, ymm8, ymm10
  vpunpcklqdq     ymm2, ymm9, ymm11  
  vpunpckhqdq     ymm3, ymm9, ymm11  

  vpmuldq      ymm4, ymm0, ymm14 // reduce(a * sub)
  vpmuldq      ymm5, ymm1, ymm14
  vpmuldq      ymm6, ymm2, ymm14
  vpmuldq      ymm7, ymm3, ymm14
  vpmuludq     ymm4, ymm4, ymm15
  vpmuludq     ymm5, ymm5, ymm15
  vpmuludq     ymm6, ymm6, ymm15
  vpmuludq     ymm7, ymm7, ymm15
  vpaddq       ymm4, ymm4, ymm0
  vpaddq       ymm5, ymm5, ymm1
  vpaddq       ymm6, ymm6, ymm2
  vpaddq       ymm7, ymm7, ymm3
  
  vpermd       ymm4, ymm12, ymm4
  vpermd       ymm5, ymm12, ymm5
  vpermd       ymm6, ymm12, ymm6
  vpermd       ymm7, ymm12, ymm7
  vmovdqa      XMMWORD PTR [reg_p1+2*rax], xmm4
  vmovdqa      XMMWORD PTR [reg_p1+2*rax+4*4], xmm5
  vmovdqa      XMMWORD PTR [reg_p1+2*rax+4*8], xmm6
  vmovdqa      XMMWORD PTR [reg_p1+2*rax+4*12], xmm7
  add          rax, 4*8
  dec          r10
  jnz          loop_pmul2
  ret

.global poly_pmul_asm3
poly_pmul_asm3:
  xor          rax, rax
  vmovdqa      ymm15, YMMWORD PTR [PARAM_Qx4+rip]
  vmovdqa      ymm14, YMMWORD PTR [PARAM_QINVx4+rip]
  vmovdqa      ymm12, YMMWORD PTR [PERM_MASK+rip]
  mov          r10, 64                        
loop_pmul3:
  vmovdqa      ymm0, [reg_p2+4*rax]                     // Load a[j]
  vmovdqa      ymm1, [reg_p2+4*rax+8*4]
  vmovdqa      ymm2, [reg_p2+4*rax+8*8]
  vmovdqa      ymm3, [reg_p2+4*rax+8*12]
  vpmovsxdq      ymm4, XMMWORD PTR [reg_p3+2*rax]                     // Load invnthroots[j]
  vpmovsxdq      ymm5, XMMWORD PTR [reg_p3+2*rax+4*4]
  vpmovsxdq      ymm6, XMMWORD PTR [reg_p3+2*rax+4*8]
  vpmovsxdq      ymm7, XMMWORD PTR [reg_p3+2*rax+4*12]
  
  vpmuldq      ymm8, ymm0, ymm4
  vpmuldq      ymm9, ymm1, ymm5
  vpmuldq      ymm10, ymm2, ymm6
  vpmuldq      ymm11, ymm3, ymm7
  vpermq       ymm4, ymm4, 0b10110001
  vpermq       ymm5, ymm5, 0b10110001
  vpermq       ymm6, ymm6, 0b10110001
  vpermq       ymm7, ymm7, 0b10110001
  vpmuldq      ymm0, ymm0, ymm4
  vpmuldq      ymm1, ymm1, ymm5
  vpmuldq      ymm2, ymm2, ymm6
  vpmuldq      ymm3, ymm3, ymm7
  
  vpunpcklqdq     ymm4, ymm8, ymm9
  vpunpckhqdq     ymm5, ymm8, ymm9
  vpunpcklqdq     ymm6, ymm10, ymm11
  vpunpckhqdq     ymm7, ymm10, ymm11
  vpsubq       ymm8, ymm4, ymm5
  vpsubq       ymm9, ymm6, ymm7

  vpunpcklqdq     ymm4, ymm0, ymm1
  vpunpckhqdq     ymm5, ymm0, ymm1
  vpunpcklqdq     ymm6, ymm2, ymm3
  vpunpckhqdq     ymm7, ymm2, ymm3
  vpaddq       ymm10, ymm4, ymm5
  vpaddq       ymm11, ymm6, ymm7
  
  vpmuldq      ymm4, ymm8, ymm14 // reduce(a * sub)
  vpmuldq      ymm5, ymm9, ymm14
  vpmuldq      ymm6, ymm10, ymm14
  vpmuldq      ymm7, ymm11, ymm14
  vpmuludq     ymm4, ymm4, ymm15
  vpmuludq     ymm5, ymm5, ymm15
  vpmuludq     ymm6, ymm6, ymm15
  vpmuludq     ymm7, ymm7, ymm15
  vpaddq       ymm4, ymm4, ymm8
  vpaddq       ymm5, ymm5, ymm9
  vpaddq       ymm6, ymm6, ymm10
  vpaddq       ymm7, ymm7, ymm11
  
  vpermq       ymm4, ymm4, 0b11011000
  vpermq       ymm5, ymm5, 0b11011000
  vpermq       ymm6, ymm6, 0b11011000
  vpermq       ymm7, ymm7, 0b11011000

  vpermd       ymm4, ymm12, ymm4
  vpermd       ymm5, ymm12, ymm5
  vpermd       ymm6, ymm12, ymm6
  vpermd       ymm7, ymm12, ymm7
  vmovdqa      XMMWORD PTR [reg_p1+rax], xmm4
  vmovdqa      XMMWORD PTR [reg_p1+rax+4*4], xmm5
  vmovdqa      XMMWORD PTR [reg_p1+rax+4*512], xmm6
  vmovdqa      XMMWORD PTR [reg_p1+rax+4*512+4*4], xmm7
  add          rax, 4*8
  dec          r10
  jnz          loop_pmul3
  ret
  
//#undef LAZYR
//**********************************************************************************
//  Inverse DGT, case N = 1024
//  Operation: c [reg_p1] <- IDGT(a) [reg_p2], 
//             [reg_p3] points to table with the inverse DGT constants.
//  Output c is in extended form. 
//********************************************************************************** 
.global poly_idgt_asm
poly_idgt_asm:
  xor          rax, rax               
  vmovdqa      ymm15, YMMWORD PTR [PARAM_Qx4+rip]
  vmovdqa      ymm14, YMMWORD PTR [PARAM_QINVx4+rip]   
  vmovdqa      ymm13, YMMWORD PTR [PARAM_BARRx4+rip]
  mov          r10, 32                               
loop_idgt1234:
// Round 1  
  vpmovsxdq    ymm0, XMMWORD PTR [reg_p2+rax]         // x[i] 
  vpmovsxdq    ymm1, XMMWORD PTR [reg_p2+rax+4*4]  
  vpmovsxdq    ymm2, XMMWORD PTR [reg_p2+rax+4*8]
  vpmovsxdq    ymm3, XMMWORD PTR [reg_p2+rax+4*12] 
  vpmovsxdq    ymm4, XMMWORD PTR [reg_p2+rax+4*16]    // x[i+m]           
  vpmovsxdq    ymm5, XMMWORD PTR [reg_p2+rax+4*20]      
  vpmovsxdq    ymm6, XMMWORD PTR [reg_p2+rax+4*24] 
  vpmovsxdq    ymm7, XMMWORD PTR [reg_p2+rax+4*28]
  
  vbroadcastss ymm12, DWORD PTR [reg_p3]
  
  vperm2i128   ymm8, ymm0, ymm1, 0x31 // sub = x[i] - x[i+m]
  vperm2i128   ymm9, ymm2, ymm3, 0x31
  vperm2i128   ymm10, ymm4, ymm5, 0x31
  vperm2i128   ymm11, ymm6, ymm7, 0x31
  vperm2i128   ymm0, ymm0, ymm1, 0x20 // sub = x[i] - x[i+m]
  vperm2i128   ymm1, ymm2, ymm3, 0x20
  vperm2i128   ymm2, ymm4, ymm5, 0x20
  vperm2i128   ymm3, ymm6, ymm7, 0x20
  vpmuldq      ymm4, ymm8, ymm12
  vpmuldq      ymm5, ymm9, ymm12
  vpmuldq      ymm6, ymm10, ymm12
  vpmuldq      ymm7, ymm11, ymm12
  
#ifdef LAZYR
  vpmuldq      ymm8, ymm4, ymm14 // reduce(a)
  vpmuldq      ymm9, ymm5, ymm14
  vpmuldq      ymm10, ymm6, ymm14
  vpmuldq      ymm11, ymm7, ymm14
  vpmuludq     ymm8, ymm8, ymm15
  vpmuludq     ymm9, ymm9, ymm15
  vpmuludq     ymm10, ymm10, ymm15
  vpmuludq     ymm11, ymm11, ymm15
  vpaddq       ymm4, ymm4, ymm8
  vpaddq       ymm5, ymm5, ymm9
  vpaddq       ymm6, ymm6, ymm10
  vpaddq       ymm7, ymm7, ymm11
  vpsrlq       ymm4, ymm4, 32
  vpsrlq       ymm5, ymm5, 32
  vpsrlq       ymm6, ymm6, 32
  vpsrlq       ymm7, ymm7, 32
#endif

  vpsubq       ymm8, ymm0, ymm4
  vpsubq       ymm9, ymm1, ymm5
  vpsubq       ymm10, ymm2, ymm6
  vpsubq       ymm11, ymm3, ymm7
  vpaddq       ymm0, ymm0, ymm4
  vpaddq       ymm1, ymm1, ymm5
  vpaddq       ymm2, ymm2, ymm6
  vpaddq       ymm3, ymm3, ymm7
  
#ifndef LAZYR
  vpmuldq      ymm4, ymm0, ymm13 // barr_reduce(add)
  vpmuldq      ymm5, ymm1, ymm13
  vpmuldq      ymm6, ymm2, ymm13
  vpmuldq      ymm7, ymm3, ymm13
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuludq     ymm4, ymm4, ymm15
  vpmuludq     ymm5, ymm5, ymm15
  vpmuludq     ymm6, ymm6, ymm15
  vpmuludq     ymm7, ymm7, ymm15
  vpsubq       ymm0, ymm0, ymm4
  vpsubq       ymm1, ymm1, ymm5
  vpsubq       ymm2, ymm2, ymm6
  vpsubq       ymm3, ymm3, ymm7

#ifndef BARRR
  vpsrad       ymm4, ymm8, 31                         // If result < 0 then add q                   
  vpsrad       ymm5, ymm9, 31                       
  vpsrad       ymm6, ymm10, 31                       
  vpsrad       ymm7, ymm11, 31   
  vpand        ymm4, ymm4, ymm15                                             
  vpand        ymm5, ymm5, ymm15                       
  vpand        ymm6, ymm6, ymm15                    
  vpand        ymm7, ymm7, ymm15
  vpaddd       ymm8, ymm8, ymm4
  vpaddd       ymm9, ymm9, ymm5
  vpaddd       ymm10, ymm10, ymm6
  vpaddd       ymm11, ymm11, ymm7
#else
  vpmuldq      ymm4, ymm8, ymm13 // barr_reduce(sub)
  vpmuldq      ymm5, ymm9, ymm13
  vpmuldq      ymm6, ymm10, ymm13
  vpmuldq      ymm7, ymm11, ymm13
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuludq     ymm4, ymm4, ymm15
  vpmuludq     ymm5, ymm5, ymm15
  vpmuludq     ymm6, ymm6, ymm15
  vpmuludq     ymm7, ymm7, ymm15
  vpsubq       ymm8, ymm8, ymm4
  vpsubq       ymm9, ymm9, ymm5
  vpsubq       ymm10, ymm10, ymm6
  vpsubq       ymm11, ymm11, ymm7
#endif
#endif

  vperm2i128   ymm4, ymm0, ymm8, 0x31 // sub = x[i] - x[i+m]
  vperm2i128   ymm5, ymm1, ymm9, 0x31
  vperm2i128   ymm6, ymm2, ymm10, 0x31
  vperm2i128   ymm7, ymm3, ymm11, 0x31
  vperm2i128   ymm0, ymm0, ymm8, 0x20 // sub = x[i] - x[i+m]
  vperm2i128   ymm1, ymm1, ymm9, 0x20
  vperm2i128   ymm2, ymm2, ymm10, 0x20
  vperm2i128   ymm3, ymm3, ymm11, 0x20
  
//Round 2
  vpmovsxdq    ymm12, XMMWORD PTR [reg_p3]
  vpmovsxdq    ymm8, XMMWORD PTR [reg_p3+4*128]
  vpunpcklqdq  ymm12, ymm12, ymm8
  vpermq       ymm12, ymm12, 0x50
  vpmuldq      ymm4, ymm4, ymm12
  vpmuldq      ymm5, ymm5, ymm12
  vpmuldq      ymm6, ymm6, ymm12
  vpmuldq      ymm7, ymm7, ymm12
  
  vpmuldq      ymm8, ymm4, ymm14 // reduce(a)
  vpmuldq      ymm9, ymm5, ymm14
  vpmuldq      ymm10, ymm6, ymm14
  vpmuldq      ymm11, ymm7, ymm14
  vpmuludq     ymm8, ymm8, ymm15
  vpmuludq     ymm9, ymm9, ymm15
  vpmuludq     ymm10, ymm10, ymm15
  vpmuludq     ymm11, ymm11, ymm15
  vpaddq       ymm4, ymm4, ymm8
  vpaddq       ymm5, ymm5, ymm9
  vpaddq       ymm6, ymm6, ymm10
  vpaddq       ymm7, ymm7, ymm11
  vpsrlq       ymm4, ymm4, 32
  vpsrlq       ymm5, ymm5, 32
  vpsrlq       ymm6, ymm6, 32
  vpsrlq       ymm7, ymm7, 32

  vpsubq       ymm8, ymm0, ymm4
  vpsubq       ymm9, ymm1, ymm5
  vpsubq       ymm10, ymm2, ymm6
  vpsubq       ymm11, ymm3, ymm7
  vpaddq       ymm0, ymm0, ymm4
  vpaddq       ymm1, ymm1, ymm5
  vpaddq       ymm2, ymm2, ymm6
  vpaddq       ymm3, ymm3, ymm7
  
#ifndef LAZYR
  vpmuldq      ymm4, ymm0, ymm13 // barr_reduce(add)
  vpmuldq      ymm5, ymm1, ymm13
  vpmuldq      ymm6, ymm2, ymm13
  vpmuldq      ymm7, ymm3, ymm13
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuludq     ymm4, ymm4, ymm15
  vpmuludq     ymm5, ymm5, ymm15
  vpmuludq     ymm6, ymm6, ymm15
  vpmuludq     ymm7, ymm7, ymm15
  vpsubq       ymm0, ymm0, ymm4
  vpsubq       ymm1, ymm1, ymm5
  vpsubq       ymm2, ymm2, ymm6
  vpsubq       ymm3, ymm3, ymm7
  
  vpsrad       ymm4, ymm8, 31                         // If result < 0 then add q                   
  vpsrad       ymm5, ymm9, 31                       
  vpsrad       ymm6, ymm10, 31                       
  vpsrad       ymm7, ymm11, 31   
  vpand        ymm4, ymm4, ymm15                                             
  vpand        ymm5, ymm5, ymm15                       
  vpand        ymm6, ymm6, ymm15                    
  vpand        ymm7, ymm7, ymm15
  vpaddd       ymm4, ymm8, ymm4
  vpaddd       ymm5, ymm9, ymm5
  vpaddd       ymm6, ymm10, ymm6
  vpaddd       ymm7, ymm11, ymm7
#else
  vmovdqa      ymm4, ymm8
  vmovdqa      ymm5, ymm9
  vmovdqa      ymm6, ymm10
  vmovdqa      ymm7, ymm11
#endif
  
//Round 3
  /*vmovdqa      [reg_p1+2*rax],     ymm0
  vmovdqa      [reg_p1+2*rax+8*4], ymm4
  vmovdqa      [reg_p1+2*rax+8*8], ymm1
  vmovdqa      [reg_p1+2*rax+8*12], ymm5
  vmovdqa      [reg_p1+2*rax+8*16], ymm2
  vmovdqa      [reg_p1+2*rax+8*20], ymm6
  vmovdqa      [reg_p1+2*rax+8*24], ymm3
  vmovdqa      [reg_p1+2*rax+8*28], ymm7*/
  
  vpmovsxdq    ymm8, XMMWORD PTR [reg_p3]
  vpmovsxdq    ymm9, XMMWORD PTR [reg_p3+4*64]
  vpmovsxdq    ymm10, XMMWORD PTR [reg_p3+4*128]
  vpmovsxdq    ymm11, XMMWORD PTR [reg_p3+4*192]
  vpunpcklqdq  ymm8, ymm8, ymm9
  vpunpcklqdq  ymm10, ymm10, ymm11
  vpermq       ymm8, ymm8, 0x50
  vpermq       ymm10, ymm10, 0x50
  vpmuldq      ymm1, ymm1, ymm8
  vpmuldq      ymm3, ymm3, ymm8
  vpmuldq      ymm5, ymm5, ymm10
  vpmuldq      ymm7, ymm7, ymm10

  vpmuldq      ymm8, ymm1, ymm14 // reduce(a)
  vpmuldq      ymm9, ymm3, ymm14
  vpmuldq      ymm10, ymm5, ymm14
  vpmuldq      ymm11, ymm7, ymm14
  vpmuludq     ymm8, ymm8, ymm15
  vpmuludq     ymm9, ymm9, ymm15
  vpmuludq     ymm10, ymm10, ymm15
  vpmuludq     ymm11, ymm11, ymm15
  vpaddq       ymm1, ymm1, ymm8
  vpaddq       ymm3, ymm3, ymm9
  vpaddq       ymm5, ymm5, ymm10
  vpaddq       ymm7, ymm7, ymm11
  vpsrlq       ymm1, ymm1, 32
  vpsrlq       ymm3, ymm3, 32
  vpsrlq       ymm5, ymm5, 32
  vpsrlq       ymm7, ymm7, 32

  vpsubq       ymm8, ymm0, ymm1
  vpsubq       ymm9, ymm4, ymm5
  vpsubq       ymm10, ymm2, ymm3
  vpsubq       ymm11, ymm6, ymm7
  vpaddq       ymm0, ymm0, ymm1
  vpaddq       ymm1, ymm4, ymm5
  vpaddq       ymm2, ymm2, ymm3
  vpaddq       ymm3, ymm6, ymm7
  
#ifndef LAZYR
  vpmuldq      ymm4, ymm0, ymm13 // barr_reduce(add)
  vpmuldq      ymm5, ymm1, ymm13
  vpmuldq      ymm6, ymm2, ymm13
  vpmuldq      ymm7, ymm3, ymm13
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuldq     ymm4, ymm4, ymm15
  vpmuldq     ymm5, ymm5, ymm15
  vpmuldq     ymm6, ymm6, ymm15
  vpmuldq     ymm7, ymm7, ymm15
  vpsubq       ymm0, ymm0, ymm4
  vpsubq       ymm1, ymm1, ymm5
  vpsubq       ymm2, ymm2, ymm6
  vpsubq       ymm3, ymm3, ymm7
  
  vpsrad       ymm4, ymm8, 31                         // If result < 0 then add q                   
  vpsrad       ymm5, ymm9, 31                       
  vpsrad       ymm6, ymm10, 31                       
  vpsrad       ymm7, ymm11, 31   
  vpand        ymm4, ymm4, ymm15                                             
  vpand        ymm5, ymm5, ymm15                       
  vpand        ymm6, ymm6, ymm15                    
  vpand        ymm7, ymm7, ymm15
  vpaddd       ymm4, ymm8, ymm4
  vpaddd       ymm5, ymm9, ymm5
  vpaddd       ymm6, ymm10, ymm6
  vpaddd       ymm7, ymm11, ymm7
#else
  vmovdqa      ymm4, ymm8
  vmovdqa      ymm5, ymm9
  vmovdqa      ymm6, ymm10
  vmovdqa      ymm7, ymm11
#endif
  
//Round 4
  /*  vmovdqa      [reg_p1+2*rax],     ymm0
  vmovdqa      [reg_p1+2*rax+8*4], ymm1
  vmovdqa      [reg_p1+2*rax+8*8], ymm4
  vmovdqa      [reg_p1+2*rax+8*12], ymm5
  vmovdqa      [reg_p1+2*rax+8*16], ymm2
  vmovdqa      [reg_p1+2*rax+8*20], ymm3
  vmovdqa      [reg_p1+2*rax+8*24], ymm6
  vmovdqa      [reg_p1+2*rax+8*28], ymm7*/
  
  vpmovsxdq    ymm8, XMMWORD PTR [reg_p3]
  vpmovsxdq    ymm9, XMMWORD PTR [reg_p3+4*32]
  vpmovsxdq    ymm10, XMMWORD PTR [reg_p3+4*64]
  vpmovsxdq    ymm11, XMMWORD PTR [reg_p3+4*96]
  vpunpcklqdq  ymm8, ymm8, ymm9
  vpunpcklqdq  ymm10, ymm10, ymm11
  vpermq       ymm8, ymm8, 0x50
  vpermq       ymm10, ymm10, 0x50
  vpmuldq      ymm2, ymm2, ymm8
  vpmuldq      ymm3, ymm3, ymm10
  vpmovsxdq    ymm8, XMMWORD PTR [reg_p3+4*128]
  vpmovsxdq    ymm9, XMMWORD PTR [reg_p3+4*160]
  vpmovsxdq    ymm10, XMMWORD PTR [reg_p3+4*192]
  vpmovsxdq    ymm11, XMMWORD PTR [reg_p3+4*224]
  vpunpcklqdq  ymm8, ymm8, ymm9
  vpunpcklqdq  ymm10, ymm10, ymm11
  vpermq       ymm8, ymm8, 0x50
  vpermq       ymm10, ymm10, 0x50
  vpmuldq      ymm6, ymm6, ymm8
  vpmuldq      ymm7, ymm7, ymm10

  vpmuldq      ymm8, ymm2, ymm14 // reduce(a)
  vpmuldq      ymm9, ymm3, ymm14
  vpmuldq      ymm10, ymm6, ymm14
  vpmuldq      ymm11, ymm7, ymm14
  vpmuludq     ymm8, ymm8, ymm15
  vpmuludq     ymm9, ymm9, ymm15
  vpmuludq     ymm10, ymm10, ymm15
  vpmuludq     ymm11, ymm11, ymm15
  vpaddq       ymm2, ymm2, ymm8
  vpaddq       ymm3, ymm3, ymm9
  vpaddq       ymm6, ymm6, ymm10
  vpaddq       ymm7, ymm7, ymm11
  vpsrlq       ymm2, ymm2, 32
  vpsrlq       ymm3, ymm3, 32
  vpsrlq       ymm6, ymm6, 32
  vpsrlq       ymm7, ymm7, 32

  vpsubq       ymm8, ymm0, ymm2
  vpsubq       ymm9, ymm1, ymm3
  vpsubq       ymm10, ymm4, ymm6
  vpsubq       ymm11, ymm5, ymm7
  vpaddq       ymm0, ymm0, ymm2
  vpaddq       ymm1, ymm1, ymm3
  vpaddq       ymm2, ymm4, ymm6
  vpaddq       ymm3, ymm5, ymm7

  vpmuldq      ymm4, ymm0, ymm13 // barr_reduce(add)
  vpmuldq      ymm5, ymm1, ymm13
  vpmuldq      ymm6, ymm2, ymm13
  vpmuldq      ymm7, ymm3, ymm13
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuldq     ymm4, ymm4, ymm15
  vpmuldq     ymm5, ymm5, ymm15
  vpmuldq     ymm6, ymm6, ymm15
  vpmuldq     ymm7, ymm7, ymm15
  vpsubq       ymm0, ymm0, ymm4
  vpsubq       ymm1, ymm1, ymm5
  vpsubq       ymm2, ymm2, ymm6
  vpsubq       ymm3, ymm3, ymm7
  
  vpmuldq      ymm4, ymm8, ymm13 // barr_reduce(sub)
  vpmuldq      ymm5, ymm9, ymm13
  vpmuldq      ymm6, ymm10, ymm13
  vpmuldq      ymm7, ymm11, ymm13
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuldq     ymm4, ymm4, ymm15
  vpmuldq     ymm5, ymm5, ymm15
  vpmuldq     ymm6, ymm6, ymm15
  vpmuldq     ymm7, ymm7, ymm15
  vpsubq       ymm4, ymm8, ymm4
  vpsubq       ymm5, ymm9, ymm5
  vpsubq       ymm6, ymm10, ymm6
  vpsubq       ymm7, ymm11, ymm7

  vmovdqa      [reg_p1+2*rax],     ymm0
  vmovdqa      [reg_p1+2*rax+8*4], ymm1
  vmovdqa      [reg_p1+2*rax+8*8], ymm2
  vmovdqa      [reg_p1+2*rax+8*12], ymm3
  vmovdqa      [reg_p1+2*rax+8*16], ymm4
  vmovdqa      [reg_p1+2*rax+8*20], ymm5
  vmovdqa      [reg_p1+2*rax+8*24], ymm6
  vmovdqa      [reg_p1+2*rax+8*28], ymm7
  
  add          rax, 4*32
  dec          r10
  jnz          loop_idgt1234

//Round 5
  xor          rax, rax
  xor          r8, r8
  xor          r11, r11
iround5_loop:
  mov          r10, 4
  vpmovsxdq    ymm8, XMMWORD PTR [reg_p3+4*r8]
  vpmovsxdq    ymm9, XMMWORD PTR [reg_p3+4*r8+4*16]
  vpunpcklqdq  ymm8, ymm8, ymm9
  vpermq       ymm12, ymm8, 0x50
loop_idgt5:
  vmovdqa      ymm0, [reg_p1+rax]
  vmovdqa      ymm1, [reg_p1+rax+8*32]
  vmovdqa      ymm2, [reg_p1+rax+8*64]
  vmovdqa      ymm3, [reg_p1+rax+8*96]
  vmovdqa      ymm4, [reg_p1+rax+8*128]
  vmovdqa      ymm5, [reg_p1+rax+8*160]
  vmovdqa      ymm6, [reg_p1+rax+8*192]
  vmovdqa      ymm7, [reg_p1+rax+8*224]

  vpmuldq      ymm1, ymm1, ymm12
  vpmuldq      ymm3, ymm3, ymm12
  vpmuldq      ymm5, ymm5, ymm12
  vpmuldq      ymm7, ymm7, ymm12

  vpmuldq      ymm8, ymm1, ymm14 // reduce(a)
  vpmuldq      ymm9, ymm3, ymm14
  vpmuldq      ymm10, ymm5, ymm14
  vpmuldq      ymm11, ymm7, ymm14
  vpmuludq     ymm8, ymm8, ymm15
  vpmuludq     ymm9, ymm9, ymm15
  vpmuludq     ymm10, ymm10, ymm15
  vpmuludq     ymm11, ymm11, ymm15
  vpaddq       ymm1, ymm1, ymm8
  vpaddq       ymm3, ymm3, ymm9
  vpaddq       ymm5, ymm5, ymm10
  vpaddq       ymm7, ymm7, ymm11
  vpsrlq       ymm1, ymm1, 32
  vpsrlq       ymm3, ymm3, 32
  vpsrlq       ymm5, ymm5, 32
  vpsrlq       ymm7, ymm7, 32

  vpsubq       ymm8, ymm0, ymm1
  vpsubq       ymm9, ymm2, ymm3
  vpsubq       ymm10, ymm4, ymm5
  vpsubq       ymm11, ymm6, ymm7
  vpaddq       ymm0, ymm0, ymm1
  vpaddq       ymm1, ymm2, ymm3
  vpaddq       ymm2, ymm4, ymm5
  vpaddq       ymm3, ymm6, ymm7

#ifndef LAZYR
  vpmuldq      ymm4, ymm0, ymm13 // barr_reduce(add)
  vpmuldq      ymm5, ymm1, ymm13
  vpmuldq      ymm6, ymm2, ymm13
  vpmuldq      ymm7, ymm3, ymm13
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuldq     ymm4, ymm4, ymm15
  vpmuldq     ymm5, ymm5, ymm15
  vpmuldq     ymm6, ymm6, ymm15
  vpmuldq     ymm7, ymm7, ymm15
  vpsubq       ymm0, ymm0, ymm4
  vpsubq       ymm1, ymm1, ymm5
  vpsubq       ymm2, ymm2, ymm6
  vpsubq       ymm3, ymm3, ymm7
  
  vpmuldq      ymm4, ymm8, ymm13 // barr_reduce(sub)
  vpmuldq      ymm5, ymm9, ymm13
  vpmuldq      ymm6, ymm10, ymm13
  vpmuldq      ymm7, ymm11, ymm13
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuldq     ymm4, ymm4, ymm15
  vpmuldq     ymm5, ymm5, ymm15
  vpmuldq     ymm6, ymm6, ymm15
  vpmuldq     ymm7, ymm7, ymm15
  vpsubq       ymm4, ymm8, ymm4
  vpsubq       ymm5, ymm9, ymm5
  vpsubq       ymm6, ymm10, ymm6
  vpsubq       ymm7, ymm11, ymm7
#else
  vmovdqa      ymm4, ymm8
  vmovdqa      ymm5, ymm9
  vmovdqa      ymm6, ymm10
  vmovdqa      ymm7, ymm11
#endif

  vmovdqa      [reg_p1+rax],     ymm0
  vmovdqa      [reg_p1+rax+8*32], ymm4
  vmovdqa      [reg_p1+rax+8*64], ymm1
  vmovdqa      [reg_p1+rax+8*96], ymm5
  vmovdqa      [reg_p1+rax+8*128], ymm2
  vmovdqa      [reg_p1+rax+8*160], ymm6
  vmovdqa      [reg_p1+rax+8*192], ymm3
  vmovdqa      [reg_p1+rax+8*224], ymm7
  
  add          rax, 8*256
  dec          r10
  jnz          loop_idgt5
  add          r8, 32
  add          r11, 4*8
  mov          rax, r11
  cmp          r11, 8*(4*8)
  jne          iround5_loop

//Round 6
  xor          rax, rax
  xor          r8, r8
  xor          r11, r11
iround6_loop:
  mov          r10, 2
  vpmovsxdq    ymm8, XMMWORD PTR [reg_p3+4*r8]
  vpmovsxdq    ymm9, XMMWORD PTR [reg_p3+4*r8+4*8]
  vpunpcklqdq  ymm8, ymm8, ymm9
  vpermq       ymm12, ymm8, 0x50
loop_idgt6:
  vmovdqa      ymm0, [reg_p1+rax]
  vmovdqa      ymm1, [reg_p1+rax+8*64]
  vmovdqa      ymm2, [reg_p1+rax+8*128]
  vmovdqa      ymm3, [reg_p1+rax+8*192]
  vmovdqa      ymm4, [reg_p1+rax+8*256]
  vmovdqa      ymm5, [reg_p1+rax+8*320]
  vmovdqa      ymm6, [reg_p1+rax+8*384]
  vmovdqa      ymm7, [reg_p1+rax+8*448]

  vpmuldq      ymm1, ymm1, ymm12
  vpmuldq      ymm3, ymm3, ymm12
  vpmuldq      ymm5, ymm5, ymm12
  vpmuldq      ymm7, ymm7, ymm12

  vpmuldq      ymm8, ymm1, ymm14 // reduce(a)
  vpmuldq      ymm9, ymm3, ymm14
  vpmuldq      ymm10, ymm5, ymm14
  vpmuldq      ymm11, ymm7, ymm14
  vpmuludq     ymm8, ymm8, ymm15
  vpmuludq     ymm9, ymm9, ymm15
  vpmuludq     ymm10, ymm10, ymm15
  vpmuludq     ymm11, ymm11, ymm15
  vpaddq       ymm1, ymm1, ymm8
  vpaddq       ymm3, ymm3, ymm9
  vpaddq       ymm5, ymm5, ymm10
  vpaddq       ymm7, ymm7, ymm11
  vpsrlq       ymm1, ymm1, 32
  vpsrlq       ymm3, ymm3, 32
  vpsrlq       ymm5, ymm5, 32
  vpsrlq       ymm7, ymm7, 32

  vpsubq       ymm8, ymm0, ymm1
  vpsubq       ymm9, ymm2, ymm3
  vpsubq       ymm10, ymm4, ymm5
  vpsubq       ymm11, ymm6, ymm7
  vpaddq       ymm0, ymm0, ymm1
  vpaddq       ymm1, ymm2, ymm3
  vpaddq       ymm2, ymm4, ymm5
  vpaddq       ymm3, ymm6, ymm7
  
#ifndef LAZYR
  vpmuldq      ymm4, ymm0, ymm13 // barr_reduce(add)
  vpmuldq      ymm5, ymm1, ymm13
  vpmuldq      ymm6, ymm2, ymm13
  vpmuldq      ymm7, ymm3, ymm13
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuldq     ymm4, ymm4, ymm15
  vpmuldq     ymm5, ymm5, ymm15
  vpmuldq     ymm6, ymm6, ymm15
  vpmuldq     ymm7, ymm7, ymm15
  vpsubq       ymm0, ymm0, ymm4
  vpsubq       ymm1, ymm1, ymm5
  vpsubq       ymm2, ymm2, ymm6
  vpsubq       ymm3, ymm3, ymm7
  
  vpmuldq      ymm4, ymm8, ymm13 // barr_reduce(sub)
  vpmuldq      ymm5, ymm9, ymm13
  vpmuldq      ymm6, ymm10, ymm13
  vpmuldq      ymm7, ymm11, ymm13
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuldq     ymm4, ymm4, ymm15
  vpmuldq     ymm5, ymm5, ymm15
  vpmuldq     ymm6, ymm6, ymm15
  vpmuldq     ymm7, ymm7, ymm15
  vpsubq       ymm4, ymm8, ymm4
  vpsubq       ymm5, ymm9, ymm5
  vpsubq       ymm6, ymm10, ymm6
  vpsubq       ymm7, ymm11, ymm7
#else
  vmovdqa      ymm4, ymm8
  vmovdqa      ymm5, ymm9
  vmovdqa      ymm6, ymm10
  vmovdqa      ymm7, ymm11
#endif

  vmovdqa      [reg_p1+rax],     ymm0
  vmovdqa      [reg_p1+rax+8*64], ymm4
  vmovdqa      [reg_p1+rax+8*128], ymm1
  vmovdqa      [reg_p1+rax+8*192], ymm5
  vmovdqa      [reg_p1+rax+8*256], ymm2
  vmovdqa      [reg_p1+rax+8*320], ymm6
  vmovdqa      [reg_p1+rax+8*384], ymm3
  vmovdqa      [reg_p1+rax+8*448], ymm7
  
  add          rax, 8*512
  dec          r10
  jnz          loop_idgt6
  add          r8, 16
  add          r11, 4*8
  mov          rax, r11
  cmp          r11, 16*(4*8)
  jne          iround6_loop

// Round 7
  xor          rax, rax 
  xor          r8, r8
  mov          r10, 32
loop_idgt789:
  vmovdqa      ymm0, [reg_p1+rax]
  vmovdqa      ymm1, [reg_p1+rax+8*128]
  vmovdqa      ymm2, [reg_p1+rax+8*256]
  vmovdqa      ymm3, [reg_p1+rax+8*384]
  vmovdqa      ymm4, [reg_p1+rax+8*512]
  vmovdqa      ymm5, [reg_p1+rax+8*640]
  vmovdqa      ymm6, [reg_p1+rax+8*768]
  vmovdqa      ymm7, [reg_p1+rax+8*896]
  
  vpmovsxdq    ymm8, XMMWORD PTR [reg_p3+4*r8]
  vpmovsxdq    ymm9, XMMWORD PTR [reg_p3+4*r8+4*4]
  vpunpcklqdq  ymm8, ymm8, ymm9
  vpermq       ymm12, ymm8, 0x50
  vpmuldq      ymm1, ymm1, ymm12
  vpmuldq      ymm3, ymm3, ymm12
  vpmuldq      ymm5, ymm5, ymm12
  vpmuldq      ymm7, ymm7, ymm12

  vpmuldq      ymm8, ymm1, ymm14 // reduce(a)
  vpmuldq      ymm9, ymm3, ymm14
  vpmuldq      ymm10, ymm5, ymm14
  vpmuldq      ymm11, ymm7, ymm14
  vpmuludq     ymm8, ymm8, ymm15
  vpmuludq     ymm9, ymm9, ymm15
  vpmuludq     ymm10, ymm10, ymm15
  vpmuludq     ymm11, ymm11, ymm15
  vpaddq       ymm1, ymm1, ymm8
  vpaddq       ymm3, ymm3, ymm9
  vpaddq       ymm5, ymm5, ymm10
  vpaddq       ymm7, ymm7, ymm11
  vpsrlq       ymm1, ymm1, 32
  vpsrlq       ymm3, ymm3, 32
  vpsrlq       ymm5, ymm5, 32
  vpsrlq       ymm7, ymm7, 32

  vpsubq       ymm8, ymm0, ymm1
  vpsubq       ymm9, ymm2, ymm3
  vpsubq       ymm10, ymm4, ymm5
  vpsubq       ymm11, ymm6, ymm7
  vpaddq       ymm0, ymm0, ymm1
  vpaddq       ymm1, ymm2, ymm3
  vpaddq       ymm2, ymm4, ymm5
  vpaddq       ymm3, ymm6, ymm7
  
#ifndef LAZYR
  vpmuldq      ymm4, ymm0, ymm13 // barr_reduce(add)
  vpmuldq      ymm5, ymm1, ymm13
  vpmuldq      ymm6, ymm2, ymm13
  vpmuldq      ymm7, ymm3, ymm13
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuldq     ymm4, ymm4, ymm15
  vpmuldq     ymm5, ymm5, ymm15
  vpmuldq     ymm6, ymm6, ymm15
  vpmuldq     ymm7, ymm7, ymm15
  vpsubq       ymm0, ymm0, ymm4
  vpsubq       ymm1, ymm1, ymm5
  vpsubq       ymm2, ymm2, ymm6
  vpsubq       ymm3, ymm3, ymm7
  
  vpmuldq      ymm4, ymm8, ymm13 // barr_reduce(sub)
  vpmuldq      ymm5, ymm9, ymm13
  vpmuldq      ymm6, ymm10, ymm13
  vpmuldq      ymm7, ymm11, ymm13
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuldq     ymm4, ymm4, ymm15
  vpmuldq     ymm5, ymm5, ymm15
  vpmuldq     ymm6, ymm6, ymm15
  vpmuldq     ymm7, ymm7, ymm15
  vpsubq       ymm4, ymm8, ymm4
  vpsubq       ymm5, ymm9, ymm5
  vpsubq       ymm6, ymm10, ymm6
  vpsubq       ymm7, ymm11, ymm7
#else
  vmovdqa      ymm4, ymm8
  vmovdqa      ymm5, ymm9
  vmovdqa      ymm6, ymm10
  vmovdqa      ymm7, ymm11
#endif

//Round 8
  vpmovsxdq    ymm12, XMMWORD PTR [reg_p3+2*r8]
  vpunpcklqdq  ymm8, ymm12, ymm12
  vpmovsxdq    ymm12, XMMWORD PTR [reg_p3+2*r8+4*128]
  vpunpcklqdq  ymm12, ymm12, ymm12
  vpmuldq      ymm1, ymm1, ymm8
  vpmuldq      ymm3, ymm3, ymm8
  vpmuldq      ymm5, ymm5, ymm12
  vpmuldq      ymm7, ymm7, ymm12

  vpmuldq      ymm8, ymm1, ymm14 // reduce(a)
  vpmuldq      ymm9, ymm3, ymm14
  vpmuldq      ymm10, ymm5, ymm14
  vpmuldq      ymm11, ymm7, ymm14
  vpmuludq     ymm8, ymm8, ymm15
  vpmuludq     ymm9, ymm9, ymm15
  vpmuludq     ymm10, ymm10, ymm15
  vpmuludq     ymm11, ymm11, ymm15
  vpaddq       ymm1, ymm1, ymm8
  vpaddq       ymm3, ymm3, ymm9
  vpaddq       ymm5, ymm5, ymm10
  vpaddq       ymm7, ymm7, ymm11
  vpsrlq       ymm1, ymm1, 32
  vpsrlq       ymm3, ymm3, 32
  vpsrlq       ymm5, ymm5, 32
  vpsrlq       ymm7, ymm7, 32

  vpsubq       ymm8, ymm0, ymm1
  vpsubq       ymm9, ymm2, ymm3
  vpsubq       ymm10, ymm4, ymm5
  vpsubq       ymm11, ymm6, ymm7
  vpaddq       ymm0, ymm0, ymm1
  vpaddq       ymm1, ymm2, ymm3
  vpaddq       ymm2, ymm4, ymm5
  vpaddq       ymm3, ymm6, ymm7
  
#ifndef LAZYR
  vpmuldq      ymm4, ymm0, ymm13 // barr_reduce(add)
  vpmuldq      ymm5, ymm1, ymm13
  vpmuldq      ymm6, ymm2, ymm13
  vpmuldq      ymm7, ymm3, ymm13
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuldq     ymm4, ymm4, ymm15
  vpmuldq     ymm5, ymm5, ymm15
  vpmuldq     ymm6, ymm6, ymm15
  vpmuldq     ymm7, ymm7, ymm15
  vpsubq       ymm0, ymm0, ymm4
  vpsubq       ymm1, ymm1, ymm5
  vpsubq       ymm2, ymm2, ymm6
  vpsubq       ymm3, ymm3, ymm7
  
  vpmuldq      ymm4, ymm8, ymm13 // barr_reduce(sub)
  vpmuldq      ymm5, ymm9, ymm13
  vpmuldq      ymm6, ymm10, ymm13
  vpmuldq      ymm7, ymm11, ymm13
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuldq     ymm4, ymm4, ymm15
  vpmuldq     ymm5, ymm5, ymm15
  vpmuldq     ymm6, ymm6, ymm15
  vpmuldq     ymm7, ymm7, ymm15
  vpsubq       ymm4, ymm8, ymm4
  vpsubq       ymm5, ymm9, ymm5
  vpsubq       ymm6, ymm10, ymm6
  vpsubq       ymm7, ymm11, ymm7
#else
  vmovdqa      ymm4, ymm8
  vmovdqa      ymm5, ymm9
  vmovdqa      ymm6, ymm10
  vmovdqa      ymm7, ymm11
#endif

// Round 9
  vpmovsxdq    ymm8, XMMWORD PTR [reg_p3+r8]
  vpmovsxdq    ymm9, XMMWORD PTR [reg_p3+r8+64*4]
  vpmovsxdq    ymm10, XMMWORD PTR [reg_p3+r8+128*4]
  vpmovsxdq    ymm11, XMMWORD PTR [reg_p3+r8+192*4]
  vpermq       ymm8, ymm8, 0x50
  vpermq       ymm9, ymm9, 0x50
  vpermq       ymm10, ymm10, 0x50
  vpermq       ymm11, ymm11, 0x50
  vpmuldq      ymm1, ymm1, ymm8
  vpmuldq      ymm3, ymm3, ymm9
  vpmuldq      ymm5, ymm5, ymm10
  vpmuldq      ymm7, ymm7, ymm11

  vpmuldq      ymm8, ymm1, ymm14 // reduce(a)
  vpmuldq      ymm9, ymm3, ymm14
  vpmuldq      ymm10, ymm5, ymm14
  vpmuldq      ymm11, ymm7, ymm14
  vpmuludq     ymm8, ymm8, ymm15
  vpmuludq     ymm9, ymm9, ymm15
  vpmuludq     ymm10, ymm10, ymm15
  vpmuludq     ymm11, ymm11, ymm15
  vpaddq       ymm1, ymm1, ymm8
  vpaddq       ymm3, ymm3, ymm9
  vpaddq       ymm5, ymm5, ymm10
  vpaddq       ymm7, ymm7, ymm11
  vpsrlq       ymm1, ymm1, 32
  vpsrlq       ymm3, ymm3, 32
  vpsrlq       ymm5, ymm5, 32
  vpsrlq       ymm7, ymm7, 32

  vpsubq       ymm8, ymm0, ymm1
  vpsubq       ymm9, ymm2, ymm3
  vpsubq       ymm10, ymm4, ymm5
  vpsubq       ymm11, ymm6, ymm7
  vpaddq       ymm0, ymm0, ymm1
  vpaddq       ymm1, ymm2, ymm3
  vpaddq       ymm2, ymm4, ymm5
  vpaddq       ymm3, ymm6, ymm7
  
  vpmuldq      ymm4, ymm0, ymm13 // barr_reduce(add)
  vpmuldq      ymm5, ymm1, ymm13
  vpmuldq      ymm6, ymm2, ymm13
  vpmuldq      ymm7, ymm3, ymm13
  vpsrlq       ymm4, ymm4, PARAM_BARR_DIV
  vpsrlq       ymm5, ymm5, PARAM_BARR_DIV
  vpsrlq       ymm6, ymm6, PARAM_BARR_DIV
  vpsrlq       ymm7, ymm7, PARAM_BARR_DIV
  vpmuldq     ymm4, ymm4, ymm15
  vpmuldq     ymm5, ymm5, ymm15
  vpmuldq     ymm6, ymm6, ymm15
  vpmuldq     ymm7, ymm7, ymm15
  vpsubq       ymm0, ymm0, ymm4
  vpsubq       ymm1, ymm1, ymm5
  vpsubq       ymm2, ymm2, ymm6
  vpsubq       ymm3, ymm3, ymm7
  
  vpsrad       ymm4, ymm8, 31                         // If result < 0 then add q                   
  vpsrad       ymm5, ymm9, 31                       
  vpsrad       ymm6, ymm10, 31                       
  vpsrad       ymm7, ymm11, 31   
  vpand        ymm4, ymm4, ymm15                                             
  vpand        ymm5, ymm5, ymm15                       
  vpand        ymm6, ymm6, ymm15                    
  vpand        ymm7, ymm7, ymm15
  vpaddd       ymm4, ymm8, ymm4
  vpaddd       ymm5, ymm9, ymm5
  vpaddd       ymm6, ymm10, ymm6
  vpaddd       ymm7, ymm11, ymm7
  
  vmovdqa      [reg_p1+rax],     ymm0
  vmovdqa      [reg_p1+rax+8*128], ymm1
  vmovdqa      [reg_p1+rax+8*256], ymm2
  vmovdqa      [reg_p1+rax+8*384], ymm3
  vmovdqa      [reg_p1+rax+8*512], ymm4
  vmovdqa      [reg_p1+rax+8*640], ymm5
  vmovdqa      [reg_p1+rax+8*768], ymm6
  vmovdqa      [reg_p1+rax+8*896], ymm7

  add          rax, 8*4
  add          r8, 8
  dec          r10
  jnz          loop_idgt789
  
  ret
