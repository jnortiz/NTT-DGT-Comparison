.intel_syntax noprefix

// Registers that are used for parameter passing:
#define reg_p1  rdi
#define reg_p2  rsi
#define reg_p3  rdx
#define reg_p4  rcx

#define PARAM_BARR_DIV 30

.text

.global poly_dgt_asm
poly_dgt_asm:

    xor          rax, rax 
    vmovdqa      ymm15, YMMWORD PTR [PARAM_Qx4+rip]
    vmovdqa      ymm14, YMMWORD PTR [PARAM_QINVx4+rip]   
    vbroadcastss ymm12, DWORD PTR [reg_p3]
    vbroadcastss ymm13, DWORD PTR [reg_p3+4]
    mov          r10, 32            
  
loop_dgt123:
    /* Round 1  */
    vmovdqa    ymm0, YMMWORD PTR [reg_p2+2*rax]         // x[i] 
    vmovdqa    ymm1, YMMWORD PTR [reg_p2+2*rax+8*128]  
    vmovdqa    ymm2, YMMWORD PTR [reg_p2+2*rax+8*256]
    vmovdqa    ymm3, YMMWORD PTR [reg_p2+2*rax+8*384] 
    vmovdqa    ymm4, YMMWORD PTR [reg_p2+2*rax+8*512]    // x[i+window]           
    vmovdqa    ymm5, YMMWORD PTR [reg_p2+2*rax+8*640]      
    vmovdqa    ymm6, YMMWORD PTR [reg_p2+2*rax+8*768] 
    vmovdqa    ymm7, YMMWORD PTR [reg_p2+2*rax+8*896]
        
    vpmuldq      ymm4, ymm4, ymm12 //a * x[i + window]
    vpmuldq      ymm5, ymm5, ymm12
    vpmuldq      ymm6, ymm6, ymm12
    vpmuldq      ymm7, ymm7, ymm12

    vpmuldq      ymm8, ymm4, ymm14 //mul = reduce(a * x[i + window])
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
    vpsrlq       ymm8, ymm4, 32                        
    vpsrlq       ymm9, ymm5, 32                       
    vpsrlq       ymm10, ymm6, 32                                      
    vpsrlq       ymm11, ymm7, 32 

    vpsubd       ymm4, ymm0, ymm8 // x[i+window] = x[i] - mul;
    vpsubd       ymm5, ymm1, ymm9
    vpsubd       ymm6, ymm2, ymm10
    vpsubd       ymm7, ymm3, ymm11
    vpaddd       ymm0, ymm0, ymm8 // x[i] = x[i] + mul
    vpaddd       ymm1, ymm1, ymm9
    vpaddd       ymm2, ymm2, ymm10
    vpaddd       ymm3, ymm3, ymm11

    vpsrad       ymm8, ymm4, 31                         // If result < 0 then add q                   
    vpsrad       ymm9, ymm5, 31                       
    vpsrad       ymm10, ymm6, 31                       
    vpsrad       ymm11, ymm7, 31   
    vpand        ymm8, ymm8, ymm15                                             
    vpand        ymm9, ymm9, ymm15                       
    vpand        ymm10, ymm10, ymm15                    
    vpand        ymm11, ymm11, ymm15
    vpaddd       ymm4, ymm4, ymm8                      
    vpaddd       ymm5, ymm5, ymm9                     
    vpaddd       ymm6, ymm6, ymm10                     
    vpaddd       ymm7, ymm7, ymm11 

    vpsubd       ymm0, ymm0, ymm15                      // If result >= q then subtract q    
    vpsubd       ymm1, ymm1, ymm15                     
    vpsubd       ymm2, ymm2, ymm15                      
    vpsubd       ymm3, ymm3, ymm15    
    vpsrad       ymm8, ymm0, 31                 
    vpsrad       ymm9, ymm1, 31                       
    vpsrad       ymm10, ymm2, 31                       
    vpsrad       ymm11, ymm3, 31   
    vpand        ymm8, ymm8, ymm15                                             
    vpand        ymm9, ymm9, ymm15                       
    vpand        ymm10, ymm10, ymm15                   
    vpand        ymm11, ymm11, ymm15
    vpaddd       ymm0, ymm0, ymm8                      
    vpaddd       ymm1, ymm1, ymm9                     
    vpaddd       ymm2, ymm2, ymm10                      
    vpaddd       ymm3, ymm3, ymm11
  
/* Round 2  */
    vpmuldq      ymm2, ymm2, ymm12 //a * x[i + window]
    vpmuldq      ymm3, ymm3, ymm12
    vpmuldq      ymm6, ymm6, ymm13
    vpmuldq      ymm7, ymm7, ymm13

    vpmuldq      ymm8, ymm2, ymm14 //mul = reduce(a * x[i + window])
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
    vpsrlq       ymm8, ymm2, 32                        
    vpsrlq       ymm9, ymm3, 32                       
    vpsrlq       ymm10, ymm6, 32                                      
    vpsrlq       ymm11, ymm7, 32 

    vpsubd       ymm2, ymm0, ymm8 // x[i+window] = x[i] - mul;
    vpsubd       ymm3, ymm1, ymm9
    vpsubd       ymm6, ymm4, ymm10
    vpsubd       ymm7, ymm5, ymm11
    vpaddd       ymm0, ymm0, ymm8 // x[i] = x[i] + mul
    vpaddd       ymm1, ymm1, ymm9
    vpaddd       ymm4, ymm4, ymm10
    vpaddd       ymm5, ymm5, ymm11

    vpsrad       ymm8, ymm2, 31                         // If result < 0 then add q                   
    vpsrad       ymm9, ymm3, 31                       
    vpsrad       ymm10, ymm6, 31                       
    vpsrad       ymm11, ymm7, 31   
    vpand        ymm8, ymm8, ymm15                                             
    vpand        ymm9, ymm9, ymm15                       
    vpand        ymm10, ymm10, ymm15                    
    vpand        ymm11, ymm11, ymm15
    vpaddd       ymm2, ymm2, ymm8                      
    vpaddd       ymm3, ymm3, ymm9                     
    vpaddd       ymm6, ymm6, ymm10                     
    vpaddd       ymm7, ymm7, ymm11 

    vpsubd       ymm0, ymm0, ymm15                      // If result >= q then subtract q    
    vpsubd       ymm1, ymm1, ymm15                     
    vpsubd       ymm4, ymm4, ymm15                      
    vpsubd       ymm5, ymm5, ymm15    
    vpsrad       ymm8, ymm0, 31                 
    vpsrad       ymm9, ymm1, 31                       
    vpsrad       ymm10, ymm4, 31                       
    vpsrad       ymm11, ymm5, 31   
    vpand        ymm8, ymm8, ymm15                                             
    vpand        ymm9, ymm9, ymm15                       
    vpand        ymm10, ymm10, ymm15                   
    vpand        ymm11, ymm11, ymm15
    vpaddd       ymm0, ymm0, ymm8                      
    vpaddd       ymm1, ymm1, ymm9                     
    vpaddd       ymm4, ymm4, ymm10                      
    vpaddd       ymm5, ymm5, ymm11
    
/* Round 3  */
    vbroadcastss ymm8, DWORD PTR [reg_p3+8]
    vbroadcastss ymm9, DWORD PTR [reg_p3+12]
    vpmuldq      ymm1, ymm1, ymm12 //a * x[i + window]
    vpmuldq      ymm3, ymm3, ymm13
    vpmuldq      ymm5, ymm5, ymm8
    vpmuldq      ymm7, ymm7, ymm9
    
    vpmuldq      ymm8, ymm1, ymm14 //mul = reduce(a * x[i + window])
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
    vpsrlq       ymm8, ymm1, 32                        
    vpsrlq       ymm9, ymm3, 32                       
    vpsrlq       ymm10, ymm5, 32                                      
    vpsrlq       ymm11, ymm7, 32 

    vpsubd       ymm1, ymm0, ymm8 // x[i+window] = x[i] - mul;
    vpsubd       ymm3, ymm2, ymm9
    vpsubd       ymm5, ymm4, ymm10
    vpsubd       ymm7, ymm6, ymm11
    vpaddd       ymm0, ymm0, ymm8 // x[i] = x[i] + mul
    vpaddd       ymm2, ymm2, ymm9
    vpaddd       ymm4, ymm4, ymm10
    vpaddd       ymm6, ymm6, ymm11

    vpsrad       ymm8, ymm1, 31                         // If result < 0 then add q                   
    vpsrad       ymm9, ymm3, 31                       
    vpsrad       ymm10, ymm5, 31                       
    vpsrad       ymm11, ymm7, 31   
    vpand        ymm8, ymm8, ymm15                                             
    vpand        ymm9, ymm9, ymm15                       
    vpand        ymm10, ymm10, ymm15                    
    vpand        ymm11, ymm11, ymm15
    vpaddd       ymm1, ymm1, ymm8                      
    vpaddd       ymm3, ymm3, ymm9                     
    vpaddd       ymm5, ymm5, ymm10                     
    vpaddd       ymm7, ymm7, ymm11 

    vpsubd       ymm0, ymm0, ymm15                      // If result >= q then subtract q    
    vpsubd       ymm2, ymm2, ymm15                     
    vpsubd       ymm4, ymm4, ymm15                      
    vpsubd       ymm6, ymm6, ymm15    
    vpsrad       ymm8, ymm0, 31                 
    vpsrad       ymm9, ymm2, 31                       
    vpsrad       ymm10, ymm4, 31                       
    vpsrad       ymm11, ymm6, 31   
    vpand        ymm8, ymm8, ymm15                                             
    vpand        ymm9, ymm9, ymm15                       
    vpand        ymm10, ymm10, ymm15                   
    vpand        ymm11, ymm11, ymm15
    vpaddd       ymm0, ymm0, ymm8                      
    vpaddd       ymm2, ymm2, ymm9                     
    vpaddd       ymm4, ymm4, ymm10                      
    vpaddd       ymm6, ymm6, ymm11
    
    vmovdqa      [reg_p1+2*rax],       ymm0
    vmovdqa      [reg_p1+2*rax+8*128], ymm1
    vmovdqa      [reg_p1+2*rax+8*256], ymm2
    vmovdqa      [reg_p1+2*rax+8*384], ymm3
    vmovdqa      [reg_p1+2*rax+8*512], ymm4
    vmovdqa      [reg_p1+2*rax+8*640], ymm5
    vmovdqa      [reg_p1+2*rax+8*768], ymm6
    vmovdqa      [reg_p1+2*rax+8*896], ymm7

    add          rax, 4*4
    dec          r10
    jnz          loop_dgt123
    
/* Round 4  */
    xor          rax, rax
    xor          r8, r8
round4loop:
    mov          r10, 4
    vbroadcastss ymm12, DWORD PTR [reg_p3+4*r8]
loop_dgt4:
    vmovdqa      ymm0, [reg_p1+rax]
    vmovdqa      ymm1, [reg_p1+rax+8*64]
    vmovdqa      ymm2, [reg_p1+rax+8*4]
    vmovdqa      ymm3, [reg_p1+rax+8*68]
    vmovdqa      ymm4, [reg_p1+rax+8*8]
    vmovdqa      ymm5, [reg_p1+rax+8*72]
    vmovdqa      ymm6, [reg_p1+rax+8*12]
    vmovdqa      ymm7, [reg_p1+rax+8*76]
    
    vpmuldq      ymm1, ymm1, ymm12
    vpmuldq      ymm3, ymm3, ymm12
    vpmuldq      ymm5, ymm5, ymm12
    vpmuldq      ymm7, ymm7, ymm12
    
    vpmuldq      ymm8, ymm1, ymm14 //mul = reduce(a * x[i + window])
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
    vpsrlq       ymm8, ymm1, 32                        
    vpsrlq       ymm9, ymm3, 32                       
    vpsrlq       ymm10, ymm5, 32                                      
    vpsrlq       ymm11, ymm7, 32 

    vpsubd       ymm1, ymm0, ymm8 // x[i+window] = x[i] - mul;
    vpsubd       ymm3, ymm2, ymm9
    vpsubd       ymm5, ymm4, ymm10
    vpsubd       ymm7, ymm6, ymm11
    vpaddd       ymm0, ymm0, ymm8 // x[i] = x[i] + mul
    vpaddd       ymm2, ymm2, ymm9
    vpaddd       ymm4, ymm4, ymm10
    vpaddd       ymm6, ymm6, ymm11

    vpsrad       ymm8, ymm1, 31                         // If result < 0 then add q                   
    vpsrad       ymm9, ymm3, 31                       
    vpsrad       ymm10, ymm5, 31                       
    vpsrad       ymm11, ymm7, 31   
    vpand        ymm8, ymm8, ymm15                                             
    vpand        ymm9, ymm9, ymm15                       
    vpand        ymm10, ymm10, ymm15                    
    vpand        ymm11, ymm11, ymm15
    vpaddd       ymm1, ymm1, ymm8                      
    vpaddd       ymm3, ymm3, ymm9                     
    vpaddd       ymm5, ymm5, ymm10                     
    vpaddd       ymm7, ymm7, ymm11 

    vpsubd       ymm0, ymm0, ymm15                      // If result >= q then subtract q    
    vpsubd       ymm2, ymm2, ymm15                     
    vpsubd       ymm4, ymm4, ymm15                      
    vpsubd       ymm6, ymm6, ymm15    
    vpsrad       ymm8, ymm0, 31                 
    vpsrad       ymm9, ymm2, 31                       
    vpsrad       ymm10, ymm4, 31                       
    vpsrad       ymm11, ymm6, 31   
    vpand        ymm8, ymm8, ymm15                                             
    vpand        ymm9, ymm9, ymm15                       
    vpand        ymm10, ymm10, ymm15                   
    vpand        ymm11, ymm11, ymm15
    vpaddd       ymm0, ymm0, ymm8                      
    vpaddd       ymm2, ymm2, ymm9                     
    vpaddd       ymm4, ymm4, ymm10                      
    vpaddd       ymm6, ymm6, ymm11

    vmovdqa      [reg_p1+rax],      ymm0
    vmovdqa      [reg_p1+rax+8*64], ymm1
    vmovdqa      [reg_p1+rax+8*4],  ymm2
    vmovdqa      [reg_p1+rax+8*68], ymm3
    vmovdqa      [reg_p1+rax+8*8],  ymm4
    vmovdqa      [reg_p1+rax+8*72], ymm5
    vmovdqa      [reg_p1+rax+8*12], ymm6
    vmovdqa      [reg_p1+rax+8*76], ymm7

    add          rax, 8*16
    dec          r10
    jnz          loop_dgt4
    add          rax, 8*64
    inc          r8
    cmp          r8, 8
    jne          round4loop

/* Round 5  */
    xor          rax, rax
    xor          r8, r8
round5loop:
    mov          r10, 2
    vbroadcastss ymm12, DWORD PTR [reg_p3+4*r8]
loop_dgt5:
    vmovdqa      ymm0, [reg_p1+rax]
    vmovdqa      ymm1, [reg_p1+rax+8*32]
    vmovdqa      ymm2, [reg_p1+rax+8*4]
    vmovdqa      ymm3, [reg_p1+rax+8*36]
    vmovdqa      ymm4, [reg_p1+rax+8*8]
    vmovdqa      ymm5, [reg_p1+rax+8*40]
    vmovdqa      ymm6, [reg_p1+rax+8*12]
    vmovdqa      ymm7, [reg_p1+rax+8*44]
    
    vpmuldq      ymm1, ymm1, ymm12
    vpmuldq      ymm3, ymm3, ymm12
    vpmuldq      ymm5, ymm5, ymm12
    vpmuldq      ymm7, ymm7, ymm12
    
    vpmuldq      ymm8, ymm1, ymm14 //mul = reduce(a * x[i + window])
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
    vpsrlq       ymm8, ymm1, 32                        
    vpsrlq       ymm9, ymm3, 32                       
    vpsrlq       ymm10, ymm5, 32                                      
    vpsrlq       ymm11, ymm7, 32 

    vpsubd       ymm1, ymm0, ymm8 // x[i+window] = x[i] - mul;
    vpsubd       ymm3, ymm2, ymm9
    vpsubd       ymm5, ymm4, ymm10
    vpsubd       ymm7, ymm6, ymm11
    vpaddd       ymm0, ymm0, ymm8 // x[i] = x[i] + mul
    vpaddd       ymm2, ymm2, ymm9
    vpaddd       ymm4, ymm4, ymm10
    vpaddd       ymm6, ymm6, ymm11

    vpsrad       ymm8, ymm1, 31                         // If result < 0 then add q                   
    vpsrad       ymm9, ymm3, 31                       
    vpsrad       ymm10, ymm5, 31                       
    vpsrad       ymm11, ymm7, 31   
    vpand        ymm8, ymm8, ymm15                                             
    vpand        ymm9, ymm9, ymm15                       
    vpand        ymm10, ymm10, ymm15                    
    vpand        ymm11, ymm11, ymm15
    vpaddd       ymm1, ymm1, ymm8                      
    vpaddd       ymm3, ymm3, ymm9                     
    vpaddd       ymm5, ymm5, ymm10                     
    vpaddd       ymm7, ymm7, ymm11 

    vpsubd       ymm0, ymm0, ymm15                      // If result >= q then subtract q    
    vpsubd       ymm2, ymm2, ymm15                     
    vpsubd       ymm4, ymm4, ymm15                      
    vpsubd       ymm6, ymm6, ymm15    
    vpsrad       ymm8, ymm0, 31                 
    vpsrad       ymm9, ymm2, 31                       
    vpsrad       ymm10, ymm4, 31                       
    vpsrad       ymm11, ymm6, 31   
    vpand        ymm8, ymm8, ymm15                                             
    vpand        ymm9, ymm9, ymm15                       
    vpand        ymm10, ymm10, ymm15                   
    vpand        ymm11, ymm11, ymm15
    vpaddd       ymm0, ymm0, ymm8                      
    vpaddd       ymm2, ymm2, ymm9                     
    vpaddd       ymm4, ymm4, ymm10                      
    vpaddd       ymm6, ymm6, ymm11

    vmovdqa      [reg_p1+rax],      ymm0
    vmovdqa      [reg_p1+rax+8*32], ymm1
    vmovdqa      [reg_p1+rax+8*4],  ymm2
    vmovdqa      [reg_p1+rax+8*36], ymm3
    vmovdqa      [reg_p1+rax+8*8],  ymm4
    vmovdqa      [reg_p1+rax+8*40], ymm5
    vmovdqa      [reg_p1+rax+8*12], ymm6
    vmovdqa      [reg_p1+rax+8*44], ymm7

    add          rax, 8*16
    dec          r10
    jnz          loop_dgt5
    add          rax, 8*32
    inc          r8
    cmp          r8, 16
    jne          round5loop

/* Round 6  */
    xor          rax, rax
    xor          r8, r8
    xor          r11, r11
loop_dgt6789:
    vmovdqa      ymm0, [reg_p1+rax]
    vmovdqa      ymm1, [reg_p1+rax+8*4]
    vmovdqa      ymm2, [reg_p1+rax+8*8]
    vmovdqa      ymm3, [reg_p1+rax+8*12]
    vmovdqa      ymm4, [reg_p1+rax+8*16]
    vmovdqa      ymm5, [reg_p1+rax+8*20]
    vmovdqa      ymm6, [reg_p1+rax+8*24]
    vmovdqa      ymm7, [reg_p1+rax+8*28]

    vbroadcastss ymm12, DWORD PTR [reg_p3+4*r8]
    vpmuldq      ymm4, ymm4, ymm12
    vpmuldq      ymm5, ymm5, ymm12
    vpmuldq      ymm6, ymm6, ymm12
    vpmuldq      ymm7, ymm7, ymm12
    
    vpmuldq      ymm8, ymm4, ymm14 //mul = reduce(a * x[i + window])
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
    vpsrlq       ymm8, ymm4, 32                        
    vpsrlq       ymm9, ymm5, 32                       
    vpsrlq       ymm10, ymm6, 32                                      
    vpsrlq       ymm11, ymm7, 32 

    vpsubd       ymm4, ymm0, ymm8 // x[i+window] = x[i] - mul;
    vpsubd       ymm5, ymm1, ymm9
    vpsubd       ymm6, ymm2, ymm10
    vpsubd       ymm7, ymm3, ymm11
    vpaddd       ymm0, ymm0, ymm8 // x[i] = x[i] + mul
    vpaddd       ymm1, ymm1, ymm9
    vpaddd       ymm2, ymm2, ymm10
    vpaddd       ymm3, ymm3, ymm11

    vpsrad       ymm8, ymm4, 31                         // If result < 0 then add q                   
    vpsrad       ymm9, ymm5, 31                       
    vpsrad       ymm10, ymm6, 31                       
    vpsrad       ymm11, ymm7, 31   
    vpand        ymm8, ymm8, ymm15                                             
    vpand        ymm9, ymm9, ymm15                       
    vpand        ymm10, ymm10, ymm15                    
    vpand        ymm11, ymm11, ymm15
    vpaddd       ymm4, ymm4, ymm8                      
    vpaddd       ymm5, ymm5, ymm9                     
    vpaddd       ymm6, ymm6, ymm10                     
    vpaddd       ymm7, ymm7, ymm11 

    vpsubd       ymm0, ymm0, ymm15                      // If result >= q then subtract q    
    vpsubd       ymm1, ymm1, ymm15                     
    vpsubd       ymm2, ymm2, ymm15                      
    vpsubd       ymm3, ymm3, ymm15    
    vpsrad       ymm8, ymm0, 31                 
    vpsrad       ymm9, ymm1, 31                       
    vpsrad       ymm10, ymm2, 31                       
    vpsrad       ymm11, ymm3, 31   
    vpand        ymm8, ymm8, ymm15                                             
    vpand        ymm9, ymm9, ymm15                       
    vpand        ymm10, ymm10, ymm15                   
    vpand        ymm11, ymm11, ymm15
    vpaddd       ymm0, ymm0, ymm8                      
    vpaddd       ymm1, ymm1, ymm9                     
    vpaddd       ymm2, ymm2, ymm10                      
    vpaddd       ymm3, ymm3, ymm11
        
/* Round 7  */
    vbroadcastss ymm8, DWORD PTR [reg_p3+8*r8]
    vbroadcastss ymm9, DWORD PTR [reg_p3+8*r8+4]
    vpmuldq      ymm2, ymm2, ymm8
    vpmuldq      ymm3, ymm3, ymm8
    vpmuldq      ymm6, ymm6, ymm9
    vpmuldq      ymm7, ymm7, ymm9
    
    vpmuldq      ymm8, ymm2, ymm14 //mul = reduce(a * x[i + window])
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
    vpsrlq       ymm8, ymm2, 32                        
    vpsrlq       ymm9, ymm3, 32                       
    vpsrlq       ymm10, ymm6, 32                                      
    vpsrlq       ymm11, ymm7, 32 

    vpsubd       ymm2, ymm0, ymm8 // x[i+window] = x[i] - mul;
    vpsubd       ymm3, ymm1, ymm9
    vpsubd       ymm6, ymm4, ymm10
    vpsubd       ymm7, ymm5, ymm11
    vpaddd       ymm0, ymm0, ymm8 // x[i] = x[i] + mul
    vpaddd       ymm1, ymm1, ymm9
    vpaddd       ymm4, ymm4, ymm10
    vpaddd       ymm5, ymm5, ymm11

    vpsrad       ymm8, ymm2, 31                         // If result < 0 then add q                   
    vpsrad       ymm9, ymm3, 31                       
    vpsrad       ymm10, ymm6, 31                       
    vpsrad       ymm11, ymm7, 31   
    vpand        ymm8, ymm8, ymm15                                             
    vpand        ymm9, ymm9, ymm15                       
    vpand        ymm10, ymm10, ymm15                    
    vpand        ymm11, ymm11, ymm15
    vpaddd       ymm2, ymm2, ymm8                      
    vpaddd       ymm3, ymm3, ymm9                     
    vpaddd       ymm6, ymm6, ymm10                     
    vpaddd       ymm7, ymm7, ymm11 

    vpsubd       ymm0, ymm0, ymm15                      // If result >= q then subtract q    
    vpsubd       ymm1, ymm1, ymm15                     
    vpsubd       ymm4, ymm4, ymm15                      
    vpsubd       ymm5, ymm5, ymm15    
    vpsrad       ymm8, ymm0, 31                 
    vpsrad       ymm9, ymm1, 31                       
    vpsrad       ymm10, ymm4, 31                       
    vpsrad       ymm11, ymm5, 31   
    vpand        ymm8, ymm8, ymm15                                             
    vpand        ymm9, ymm9, ymm15                       
    vpand        ymm10, ymm10, ymm15                   
    vpand        ymm11, ymm11, ymm15
    vpaddd       ymm0, ymm0, ymm8                      
    vpaddd       ymm1, ymm1, ymm9                     
    vpaddd       ymm4, ymm4, ymm10                      
    vpaddd       ymm5, ymm5, ymm11

/* Round 8  */
    vbroadcastss ymm8, DWORD PTR [reg_p3+4*r11]
    vbroadcastss ymm9, DWORD PTR [reg_p3+4*r11+4]
    vbroadcastss ymm10, DWORD PTR [reg_p3+4*r11+8]
    vbroadcastss ymm11, DWORD PTR [reg_p3+4*r11+12]
    vpmuldq      ymm1, ymm1, ymm8
    vpmuldq      ymm3, ymm3, ymm9
    vpmuldq      ymm5, ymm5, ymm10
    vpmuldq      ymm7, ymm7, ymm11
    
    vpmuldq      ymm8, ymm1, ymm14 //mul = reduce(a * x[i + window])
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
    vpsrlq       ymm8, ymm1, 32                        
    vpsrlq       ymm9, ymm3, 32                       
    vpsrlq       ymm10, ymm5, 32                                      
    vpsrlq       ymm11, ymm7, 32 

    vpsubd       ymm1, ymm0, ymm8 // x[i+window] = x[i] - mul;
    vpsubd       ymm3, ymm2, ymm9
    vpsubd       ymm5, ymm4, ymm10
    vpsubd       ymm7, ymm6, ymm11
    vpaddd       ymm0, ymm0, ymm8 // x[i] = x[i] + mul
    vpaddd       ymm2, ymm2, ymm9
    vpaddd       ymm4, ymm4, ymm10
    vpaddd       ymm6, ymm6, ymm11

    vpsrad       ymm8, ymm1, 31                         // If result < 0 then add q                   
    vpsrad       ymm9, ymm3, 31                       
    vpsrad       ymm10, ymm5, 31                       
    vpsrad       ymm11, ymm7, 31   
    vpand        ymm8, ymm8, ymm15                                             
    vpand        ymm9, ymm9, ymm15                       
    vpand        ymm10, ymm10, ymm15                    
    vpand        ymm11, ymm11, ymm15
    vpaddd       ymm1, ymm1, ymm8                      
    vpaddd       ymm3, ymm3, ymm9                     
    vpaddd       ymm5, ymm5, ymm10                     
    vpaddd       ymm7, ymm7, ymm11 

    vpsubd       ymm0, ymm0, ymm15                      // If result >= q then subtract q    
    vpsubd       ymm2, ymm2, ymm15                     
    vpsubd       ymm4, ymm4, ymm15                      
    vpsubd       ymm6, ymm6, ymm15    
    vpsrad       ymm8, ymm0, 31                 
    vpsrad       ymm9, ymm2, 31                       
    vpsrad       ymm10, ymm4, 31                       
    vpsrad       ymm11, ymm6, 31   
    vpand        ymm8, ymm8, ymm15                                             
    vpand        ymm9, ymm9, ymm15                       
    vpand        ymm10, ymm10, ymm15                   
    vpand        ymm11, ymm11, ymm15
    vpaddd       ymm0, ymm0, ymm8                      
    vpaddd       ymm2, ymm2, ymm9                     
    vpaddd       ymm4, ymm4, ymm10                      
    vpaddd       ymm6, ymm6, ymm11
    
/* Round 9  */
    vperm2i128   ymm12, ymm0, ymm1, 0x20
    vperm2i128   ymm0, ymm0, ymm1, 0x31
    vperm2i128   ymm1, ymm2, ymm3, 0x20
    vperm2i128   ymm2, ymm2, ymm3, 0x31
    vperm2i128   ymm3, ymm4, ymm5, 0x20
    vperm2i128   ymm4, ymm4, ymm5, 0x31
    vperm2i128   ymm5, ymm6, ymm7, 0x20
    vperm2i128   ymm6, ymm6, ymm7, 0x31
    
    vpmovzxdq    ymm8, XMMWORD PTR [reg_p3+8*r11]
    vpmovzxdq    ymm9, XMMWORD PTR [reg_p3+8*r11+16]
    vpermq       ymm10, ymm8, 0x50
    vpermq       ymm11, ymm9, 0x50
    vpmuldq      ymm0, ymm0, ymm10
    vpmuldq      ymm4, ymm4, ymm11
    vpermq       ymm10, ymm8, 0xFA
    vpermq       ymm11, ymm9, 0xFA
    vpmuldq      ymm2, ymm2, ymm10
    vpmuldq      ymm6, ymm6, ymm11
    
    vpmuldq      ymm8, ymm0, ymm14 //mul = reduce(a * x[i + window])
    vpmuldq      ymm9, ymm2, ymm14                    
    vpmuldq      ymm10, ymm4, ymm14 
    vpmuldq      ymm11, ymm6, ymm14
    vpmuludq     ymm8, ymm8, ymm15  
    vpmuludq     ymm9, ymm9, ymm15 
    vpmuludq     ymm10, ymm10, ymm15                   
    vpmuludq     ymm11, ymm11, ymm15
    vpaddq       ymm0, ymm0, ymm8                      
    vpaddq       ymm2, ymm2, ymm9                     
    vpaddq       ymm4, ymm4, ymm10                                               
    vpaddq       ymm6, ymm6, ymm11  
    vpsrlq       ymm8, ymm0, 32                        
    vpsrlq       ymm9, ymm2, 32                       
    vpsrlq       ymm10, ymm4, 32                                      
    vpsrlq       ymm11, ymm6, 32 

    vpsubd       ymm0, ymm12, ymm8 // x[i+window] = x[i] - mul;
    vpsubd       ymm2, ymm1, ymm9
    vpsubd       ymm4, ymm3, ymm10
    vpsubd       ymm6, ymm5, ymm11
    vpaddd       ymm12, ymm12, ymm8 // x[i] = x[i] + mul
    vpaddd       ymm1, ymm1, ymm9
    vpaddd       ymm3, ymm3, ymm10
    vpaddd       ymm5, ymm5, ymm11

    vpsrad       ymm8, ymm0, 31                         // If result < 0 then add q                   
    vpsrad       ymm9, ymm2, 31                       
    vpsrad       ymm10, ymm4, 31                       
    vpsrad       ymm11, ymm6, 31   
    vpand        ymm8, ymm8, ymm15                                             
    vpand        ymm9, ymm9, ymm15                       
    vpand        ymm10, ymm10, ymm15                    
    vpand        ymm11, ymm11, ymm15
    vpaddd       ymm0, ymm0, ymm8                      
    vpaddd       ymm2, ymm2, ymm9                     
    vpaddd       ymm4, ymm4, ymm10                     
    vpaddd       ymm6, ymm6, ymm11 

    vpsubd       ymm12, ymm12, ymm15                      // If result >= q then subtract q    
    vpsubd       ymm1, ymm1, ymm15                     
    vpsubd       ymm3, ymm3, ymm15                      
    vpsubd       ymm5, ymm5, ymm15    
    vpsrad       ymm8, ymm12, 31                 
    vpsrad       ymm9, ymm1, 31                       
    vpsrad       ymm10, ymm3, 31                       
    vpsrad       ymm11, ymm5, 31   
    vpand        ymm8, ymm8, ymm15                                             
    vpand        ymm9, ymm9, ymm15                       
    vpand        ymm10, ymm10, ymm15                   
    vpand        ymm11, ymm11, ymm15
    vpaddd       ymm12, ymm12, ymm8                      
    vpaddd       ymm1, ymm1, ymm9                     
    vpaddd       ymm3, ymm3, ymm10                      
    vpaddd       ymm5, ymm5, ymm11
    
    vperm2i128   ymm7, ymm12, ymm0, 0x20
    vperm2i128   ymm12, ymm12, ymm0, 0x31
    vperm2i128   ymm0, ymm1, ymm2, 0x20
    vperm2i128   ymm1, ymm1, ymm2, 0x31
    vperm2i128   ymm2, ymm3, ymm4, 0x20
    vperm2i128   ymm3, ymm3, ymm4, 0x31
    vperm2i128   ymm4, ymm5, ymm6, 0x20
    vperm2i128   ymm5, ymm5, ymm6, 0x31
    
    vmovdqa      [reg_p1+rax],      ymm7
    vmovdqa      [reg_p1+rax+8*4],  ymm12
    vmovdqa      [reg_p1+rax+8*8],  ymm0
    vmovdqa      [reg_p1+rax+8*12], ymm1
    vmovdqa      [reg_p1+rax+8*16], ymm2
    vmovdqa      [reg_p1+rax+8*20], ymm3
    vmovdqa      [reg_p1+rax+8*24], ymm4
    vmovdqa      [reg_p1+rax+8*28], ymm5

    add          rax, 8*32
    inc          r8
    add          r11, 4
    cmp          r8, 32
    jne          loop_dgt6789
    
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
  mov          r10, 64                               
loop_pmul:
  vmovdqa      ymm0, [reg_p2+rax]                     // Load a[j]
  vmovdqa      ymm2, [reg_p2+rax+4*512]
  vmovdqa      ymm4, [reg_p3+rax]                      
  vmovdqa      ymm6, [reg_p3+rax+4*512]
  
  vpsrlq       ymm1, ymm0, 32
  vpsrlq       ymm3, ymm2, 32
  vpsrlq       ymm5, ymm4, 32
  vpsrlq       ymm7, ymm6, 32
  
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
  
  vperm2i128   ymm0, ymm8, ymm9, 0x20
  vperm2i128   ymm1, ymm8, ymm9, 0x31
  vperm2i128   ymm2, ymm10, ymm11, 0x20
  vperm2i128   ymm3, ymm10, ymm11, 0x31

  vpunpcklqdq     ymm8, ymm0, ymm2
  vpunpcklqdq     ymm9, ymm1, ymm3
  vpunpckhqdq     ymm10, ymm0, ymm2
  vpunpckhqdq     ymm11, ymm1, ymm3  
  vpmuldq      ymm4, ymm8, ymm14
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
  vpsrlq       ymm4, ymm4, 32
  vpsrlq       ymm5, ymm5, 32
  vpsrlq       ymm6, ymm6, 32
  vpsrlq       ymm7, ymm7, 32
  vmovdqa      [reg_p1+4*rax], ymm4
  vmovdqa      [reg_p1+4*rax+8*4], ymm6
  vmovdqa      [reg_p1+4*rax+8*8], ymm5
  vmovdqa      [reg_p1+4*rax+8*12], ymm7

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
  mov          r10, 64                              
loop_pmul2:
  vmovdqa      ymm0, YMMWORD PTR [reg_p2+rax]                     // Load a[j]
  vmovdqa      ymm2, YMMWORD PTR [reg_p2+rax+4*512]
  vmovdqa      ymm4, YMMWORD PTR [reg_p3+4*rax]
  vmovdqa      ymm5, YMMWORD PTR [reg_p3+4*rax+8*4]
  vmovdqa      ymm6, YMMWORD PTR [reg_p3+4*rax+8*8]
  vmovdqa      ymm7, YMMWORD PTR [reg_p3+4*rax+8*12]
  
  vpsrlq        ymm1, ymm0, 32
  vpsrlq        ymm3, ymm2, 32
  
  vpunpcklqdq     ymm8, ymm4, ymm5
  vpunpckhqdq     ymm9, ymm4, ymm5
  vpunpcklqdq     ymm10, ymm6, ymm7
  vpunpckhqdq     ymm11, ymm6, ymm7
  
  vpermq       ymm0, ymm0, 0b11011000
  vpermq       ymm1, ymm1, 0b11011000
  vpermq       ymm2, ymm2, 0b11011000
  vpermq       ymm3, ymm3, 0b11011000
  
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

  vpunpcklqdq     ymm0, ymm8, ymm10
  vpunpckhqdq     ymm1, ymm8, ymm10
  vpunpcklqdq     ymm2, ymm9, ymm11  
  vpunpckhqdq     ymm3, ymm9, ymm11  

  vpmuldq      ymm4, ymm0, ymm14 // reduce
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
  vpsrlq       ymm4, ymm4, 32
  vpsrlq       ymm5, ymm5, 32
  vpsrlq       ymm6, ymm6, 32
  vpsrlq       ymm7, ymm7, 32
  
  vmovdqa      [reg_p1+4*rax], ymm4
  vmovdqa      [reg_p1+4*rax+8*4], ymm5
  vmovdqa      [reg_p1+4*rax+8*8], ymm6
  vmovdqa      [reg_p1+4*rax+8*12], ymm7
  add          rax, 4*8
  dec          r10
  jnz          loop_pmul2
  ret

    
.global poly_pmul_asm3
poly_pmul_asm3:
  xor          rax, rax
  vmovdqa      ymm15, YMMWORD PTR [PARAM_Qx4+rip]
  vmovdqa      ymm14, YMMWORD PTR [PARAM_QINVx4+rip]
  vmovdqa      ymm12, YMMWORD PTR [EXIT_MASK+rip]
  mov          r10, 64                        
loop_pmul3:
  vmovdqa      ymm0, [reg_p2+4*rax]                     // Load a[j]
  vmovdqa      ymm1, [reg_p2+4*rax+8*4]
  vmovdqa      ymm2, [reg_p2+4*rax+8*8]
  vmovdqa      ymm3, [reg_p2+4*rax+8*12]
  vmovdqa      ymm4, [reg_p3+2*rax]                     // Load invnthroots[j]
  vmovdqa      ymm6, [reg_p3+2*rax+8*4]
  vpsrlq       ymm5, ymm4, 32
  vpsrlq       ymm7, ymm6, 32
  
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
  
  vpermd       ymm4, ymm12, ymm4
  vpermd       ymm5, ymm12, ymm5
  vpermd       ymm6, ymm12, ymm6
  vpermd       ymm7, ymm12, ymm7
  vperm2i128   ymm4, ymm4, ymm5, 0x20
  vperm2i128   ymm6, ymm6, ymm7, 0x20
  vmovdqa      [reg_p1+rax], ymm4
  vmovdqa      [reg_p1+rax+4*512], ymm6

  add          rax, 4*8
  dec          r10
  jnz          loop_pmul3
  ret

//**********************************************************************************
//  Inverse DGT, case N = 1024
//  Operation: c [reg_p1] <- IDGT(a) [reg_p2], 
//             [reg_p3] points to table with the inverse DGT constants.
//  Output c is in extended form. 
//********************************************************************************** 
.global poly_idgt_asm
poly_idgt_asm:
    xor         rax, rax               
    vmovdqa     ymm15, YMMWORD PTR [PARAM_Qx4+rip]
    vmovdqa     ymm14, YMMWORD PTR [PARAM_QINVx4+rip]   
    vmovdqa     ymm13, YMMWORD PTR [PARAM_BARRx4+rip]
    mov         r10, 32
    xor         r11, r11
    xor         r8, r8
  
loop_idgt1234:
/* Round 1  */
    vmovdqa     ymm0, [reg_p2+2*rax]      // x[i]
    vmovdqa     ymm1, [reg_p2+2*rax+8*4]
    vmovdqa     ymm2, [reg_p2+2*rax+8*8]
    vmovdqa     ymm3, [reg_p2+2*rax+8*12]
    vmovdqa     ymm4, [reg_p2+2*rax+8*16] // x[i + window]
    vmovdqa     ymm5, [reg_p2+2*rax+8*20]
    vmovdqa     ymm6, [reg_p2+2*rax+8*24]
    vmovdqa     ymm7, [reg_p2+2*rax+8*28]
    
    vperm2i128  ymm8, ymm0, ymm1, 0x31
    vperm2i128  ymm0, ymm0, ymm1, 0x20
    vperm2i128  ymm9, ymm2, ymm3, 0x31
    vperm2i128  ymm1, ymm2, ymm3, 0x20
    vperm2i128  ymm10, ymm4, ymm5, 0x31
    vperm2i128  ymm2, ymm4, ymm5, 0x20
    vperm2i128  ymm11, ymm6, ymm7, 0x31
    vperm2i128  ymm3, ymm6, ymm7, 0x20
    
    vpsubd      ymm4, ymm0, ymm8
    vpsubd      ymm5, ymm1, ymm9
    vpsubd      ymm6, ymm2, ymm10
    vpsubd      ymm7, ymm3, ymm11
    
    vpaddd      ymm0, ymm0, ymm8
    vpaddd      ymm1, ymm1, ymm9
    vpaddd      ymm2, ymm2, ymm10
    vpaddd      ymm3, ymm3, ymm11

    vpmuldq     ymm8, ymm0, ymm13 // barr_reduce(add)
    vpmuldq     ymm9, ymm1, ymm13
    vpmuldq     ymm10, ymm2, ymm13
    vpmuldq     ymm11, ymm3, ymm13
    vpsrlq      ymm8, ymm8, PARAM_BARR_DIV
    vpsrlq      ymm9, ymm9, PARAM_BARR_DIV
    vpsrlq      ymm10, ymm10, PARAM_BARR_DIV
    vpsrlq      ymm11, ymm11, PARAM_BARR_DIV
    vpmuludq    ymm8, ymm8, ymm15
    vpmuludq    ymm9, ymm9, ymm15
    vpmuludq    ymm10, ymm10, ymm15
    vpmuludq    ymm11, ymm11, ymm15
    vpsubq      ymm0, ymm0, ymm8
    vpsubq      ymm1, ymm1, ymm9
    vpsubq      ymm2, ymm2, ymm10
    vpsubq      ymm3, ymm3, ymm11
    
    vpmovzxdq   ymm12, XMMWORD PTR [reg_p3+8*r11]
    vpermq      ymm8, ymm12, 0b01010000
    vpmuldq     ymm4, ymm4, ymm8
    vpermq      ymm9, ymm12, 0b11111010
    vpmuldq     ymm5, ymm5, ymm9
    vpmovzxdq   ymm12, XMMWORD PTR [reg_p3+8*r11 + 4*4]
    vpermq      ymm10, ymm12, 0b01010000
    vpmuldq     ymm6, ymm6, ymm10
    vpermq      ymm11, ymm12, 0b11111010
    vpmuldq     ymm7, ymm7, ymm11
    
    vpmuldq     ymm8, ymm4, ymm14 //mul = reduce(a * sub)
    vpmuldq     ymm9, ymm5, ymm14                    
    vpmuldq     ymm10, ymm6, ymm14 
    vpmuldq     ymm11, ymm7, ymm14
    vpmuludq    ymm8, ymm8, ymm15  
    vpmuludq    ymm9, ymm9, ymm15 
    vpmuludq    ymm10, ymm10, ymm15                   
    vpmuludq    ymm11, ymm11, ymm15
    vpaddq      ymm4, ymm4, ymm8                      
    vpaddq      ymm5, ymm5, ymm9                     
    vpaddq      ymm6, ymm6, ymm10                                               
    vpaddq      ymm7, ymm7, ymm11  
    vpsrlq      ymm8, ymm4, 32                        
    vpsrlq      ymm9, ymm5, 32                       
    vpsrlq      ymm10, ymm6, 32                                      
    vpsrlq      ymm11, ymm7, 32 

    vperm2i128  ymm4, ymm0, ymm8, 0x31
    vperm2i128  ymm5, ymm1, ymm9, 0x31
    vperm2i128  ymm6, ymm2, ymm10, 0x31
    vperm2i128  ymm7, ymm3, ymm11, 0x31
    vperm2i128  ymm0, ymm0, ymm8, 0x20
    vperm2i128  ymm1, ymm1, ymm9, 0x20
    vperm2i128  ymm2, ymm2, ymm10, 0x20
    vperm2i128  ymm3, ymm3, ymm11, 0x20
    
    vpsubd      ymm8, ymm0, ymm4
    vpsubd      ymm9, ymm1, ymm5
    vpsubd      ymm10, ymm2, ymm6
    vpsubd      ymm11, ymm3, ymm7
    
    vpaddd      ymm0, ymm0, ymm4
    vpaddd      ymm1, ymm1, ymm5
    vpaddd      ymm2, ymm2, ymm6
    vpaddd      ymm3, ymm3, ymm7

    vpmuldq     ymm4, ymm0, ymm13 // barr_reduce(add)
    vpmuldq     ymm5, ymm1, ymm13
    vpmuldq     ymm6, ymm2, ymm13
    vpmuldq     ymm7, ymm3, ymm13
    vpsrlq      ymm4, ymm4, PARAM_BARR_DIV
    vpsrlq      ymm5, ymm5, PARAM_BARR_DIV
    vpsrlq      ymm6, ymm6, PARAM_BARR_DIV
    vpsrlq      ymm7, ymm7, PARAM_BARR_DIV
    vpmuludq    ymm4, ymm4, ymm15
    vpmuludq    ymm5, ymm5, ymm15
    vpmuludq    ymm6, ymm6, ymm15
    vpmuludq    ymm7, ymm7, ymm15
    vpsubq      ymm0, ymm0, ymm4
    vpsubq      ymm1, ymm1, ymm5
    vpsubq      ymm2, ymm2, ymm6
    vpsubq      ymm3, ymm3, ymm7
    
    vbroadcastss ymm12, DWORD PTR [reg_p3 + 4*r11]
    vpmuldq     ymm4, ymm8, ymm12
    vbroadcastss ymm12, DWORD PTR [reg_p3 + 4*r11 + 4]
    vpmuldq     ymm5, ymm9, ymm12
    vbroadcastss ymm12, DWORD PTR [reg_p3 + 4*r11 + 8]
    vpmuldq     ymm6, ymm10, ymm12
    vbroadcastss ymm12, DWORD PTR [reg_p3 + 4*r11 + 12]
    vpmuldq     ymm7, ymm11, ymm12
    
    vpmuldq     ymm8, ymm4, ymm14 //mul = reduce(a * sub)
    vpmuldq     ymm9, ymm5, ymm14                    
    vpmuldq     ymm10, ymm6, ymm14 
    vpmuldq     ymm11, ymm7, ymm14
    vpmuludq    ymm8, ymm8, ymm15  
    vpmuludq    ymm9, ymm9, ymm15 
    vpmuludq    ymm10, ymm10, ymm15                   
    vpmuludq    ymm11, ymm11, ymm15
    vpaddq      ymm4, ymm4, ymm8                      
    vpaddq      ymm5, ymm5, ymm9                     
    vpaddq      ymm6, ymm6, ymm10                                               
    vpaddq      ymm7, ymm7, ymm11  
    vpsrlq      ymm4, ymm4, 32                        
    vpsrlq      ymm5, ymm5, 32                       
    vpsrlq      ymm6, ymm6, 32                                      
    vpsrlq      ymm7, ymm7, 32 
    
/* Round 3 */
    vpsubd      ymm8, ymm0, ymm1
    vpsubd      ymm9, ymm2, ymm3
    vpsubd      ymm10, ymm4, ymm5
    vpsubd      ymm11, ymm6, ymm7
    
    vpaddd      ymm0, ymm0, ymm1
    vpaddd      ymm1, ymm2, ymm3
    vpaddd      ymm2, ymm4, ymm5
    vpaddd      ymm3, ymm6, ymm7

    vpmuldq     ymm4, ymm0, ymm13 // barr_reduce(add)
    vpmuldq     ymm5, ymm1, ymm13
    vpmuldq     ymm6, ymm2, ymm13
    vpmuldq     ymm7, ymm3, ymm13
    vpsrlq      ymm4, ymm4, PARAM_BARR_DIV
    vpsrlq      ymm5, ymm5, PARAM_BARR_DIV
    vpsrlq      ymm6, ymm6, PARAM_BARR_DIV
    vpsrlq      ymm7, ymm7, PARAM_BARR_DIV
    vpmuludq    ymm4, ymm4, ymm15
    vpmuludq    ymm5, ymm5, ymm15
    vpmuludq    ymm6, ymm6, ymm15
    vpmuludq    ymm7, ymm7, ymm15
    vpsubq      ymm0, ymm0, ymm4
    vpsubq      ymm1, ymm1, ymm5
    vpsubq      ymm2, ymm2, ymm6
    vpsubq      ymm3, ymm3, ymm7
    
    vbroadcastss ymm12, DWORD PTR [reg_p3+8*r8]
    vpmuldq     ymm4, ymm8, ymm12
    vpmuldq     ymm6, ymm10, ymm12
    vbroadcastss ymm12, DWORD PTR [reg_p3+8*r8+4]
    vpmuldq     ymm5, ymm9, ymm12
    vpmuldq     ymm7, ymm11, ymm12
    
    vpmuldq     ymm8, ymm4, ymm14 //mul = reduce(a * sub)
    vpmuldq     ymm9, ymm5, ymm14                    
    vpmuldq     ymm10, ymm6, ymm14 
    vpmuldq     ymm11, ymm7, ymm14
    vpmuludq    ymm8, ymm8, ymm15  
    vpmuludq    ymm9, ymm9, ymm15 
    vpmuludq    ymm10, ymm10, ymm15                   
    vpmuludq    ymm11, ymm11, ymm15
    vpaddq      ymm4, ymm4, ymm8                      
    vpaddq      ymm5, ymm5, ymm9                     
    vpaddq      ymm6, ymm6, ymm10                                               
    vpaddq      ymm7, ymm7, ymm11  
    vpsrlq      ymm4, ymm4, 32                        
    vpsrlq      ymm5, ymm5, 32                       
    vpsrlq      ymm6, ymm6, 32                                      
    vpsrlq      ymm7, ymm7, 32 
    
/* Round 4 */
    vpsubd      ymm8, ymm0, ymm1
    vpsubd      ymm9, ymm2, ymm3
    vpsubd      ymm10, ymm4, ymm5
    vpsubd      ymm11, ymm6, ymm7
    
    vpaddd      ymm0, ymm0, ymm1
    vpaddd      ymm1, ymm2, ymm3
    vpaddd      ymm2, ymm4, ymm5
    vpaddd      ymm3, ymm6, ymm7

    vpmuldq     ymm4, ymm0, ymm13 // barr_reduce(add)
    vpmuldq     ymm5, ymm1, ymm13
    vpmuldq     ymm6, ymm2, ymm13
    vpmuldq     ymm7, ymm3, ymm13
    vpsrlq      ymm4, ymm4, PARAM_BARR_DIV
    vpsrlq      ymm5, ymm5, PARAM_BARR_DIV
    vpsrlq      ymm6, ymm6, PARAM_BARR_DIV
    vpsrlq      ymm7, ymm7, PARAM_BARR_DIV
    vpmuludq    ymm4, ymm4, ymm15
    vpmuludq    ymm5, ymm5, ymm15
    vpmuludq    ymm6, ymm6, ymm15
    vpmuludq    ymm7, ymm7, ymm15
    vpsubq      ymm0, ymm0, ymm4
    vpsubq      ymm1, ymm1, ymm5
    vpsubq      ymm2, ymm2, ymm6
    vpsubq      ymm3, ymm3, ymm7
    
    vbroadcastss ymm12, DWORD PTR [reg_p3+4*r8]
    vpmuldq     ymm4, ymm8, ymm12
    vpmuldq     ymm5, ymm9, ymm12
    vpmuldq     ymm6, ymm10, ymm12
    vpmuldq     ymm7, ymm11, ymm12
    
    vpmuldq     ymm8, ymm4, ymm14 //mul = reduce(a * sub)
    vpmuldq     ymm9, ymm5, ymm14                    
    vpmuldq     ymm10, ymm6, ymm14 
    vpmuldq     ymm11, ymm7, ymm14
    vpmuludq    ymm8, ymm8, ymm15  
    vpmuludq    ymm9, ymm9, ymm15 
    vpmuludq    ymm10, ymm10, ymm15                   
    vpmuludq    ymm11, ymm11, ymm15
    vpaddq      ymm4, ymm4, ymm8                      
    vpaddq      ymm5, ymm5, ymm9                     
    vpaddq      ymm6, ymm6, ymm10                                               
    vpaddq      ymm7, ymm7, ymm11  
    vpsrlq      ymm4, ymm4, 32                        
    vpsrlq      ymm5, ymm5, 32                       
    vpsrlq      ymm6, ymm6, 32                                      
    vpsrlq      ymm7, ymm7, 32
        
    vmovdqa     [reg_p1+2*rax],      ymm0
    vmovdqa     [reg_p1+2*rax+8*4],  ymm1
    vmovdqa     [reg_p1+2*rax+8*8],  ymm2
    vmovdqa     [reg_p1+2*rax+8*12], ymm3
    vmovdqa     [reg_p1+2*rax+8*16], ymm4
    vmovdqa     [reg_p1+2*rax+8*20], ymm5
    vmovdqa     [reg_p1+2*rax+8*24], ymm6
    vmovdqa     [reg_p1+2*rax+8*28], ymm7
    
    add         rax, 4*32
    add         r11, 4
    inc         r8
    dec         r10
    jnz         loop_idgt1234
    
/* Round 5  */
    xor          rax, rax
    xor          r8, r8
iround5loop:
    mov          r10, 2
    vbroadcastss ymm12, DWORD PTR [reg_p3+4*r8]
loop_idgt5:
    vmovdqa      ymm0, [reg_p1+rax]
    vmovdqa      ymm1, [reg_p1+rax+8*32]
    vmovdqa      ymm2, [reg_p1+rax+8*4]
    vmovdqa      ymm3, [reg_p1+rax+8*36]
    vmovdqa      ymm4, [reg_p1+rax+8*8]
    vmovdqa      ymm5, [reg_p1+rax+8*40]
    vmovdqa      ymm6, [reg_p1+rax+8*12]
    vmovdqa      ymm7, [reg_p1+rax+8*44]
    
    vpsubd      ymm8, ymm0, ymm1
    vpsubd      ymm9, ymm2, ymm3
    vpsubd      ymm10, ymm4, ymm5
    vpsubd      ymm11, ymm6, ymm7
    
    vpaddd      ymm0, ymm0, ymm1
    vpaddd      ymm1, ymm2, ymm3
    vpaddd      ymm2, ymm4, ymm5
    vpaddd      ymm3, ymm6, ymm7

    vpmuldq     ymm4, ymm0, ymm13 // barr_reduce(add)
    vpmuldq     ymm5, ymm1, ymm13
    vpmuldq     ymm6, ymm2, ymm13
    vpmuldq     ymm7, ymm3, ymm13
    vpsrlq      ymm4, ymm4, PARAM_BARR_DIV
    vpsrlq      ymm5, ymm5, PARAM_BARR_DIV
    vpsrlq      ymm6, ymm6, PARAM_BARR_DIV
    vpsrlq      ymm7, ymm7, PARAM_BARR_DIV
    vpmuludq    ymm4, ymm4, ymm15
    vpmuludq    ymm5, ymm5, ymm15
    vpmuludq    ymm6, ymm6, ymm15
    vpmuludq    ymm7, ymm7, ymm15
    vpsubq      ymm0, ymm0, ymm4
    vpsubq      ymm1, ymm1, ymm5
    vpsubq      ymm2, ymm2, ymm6
    vpsubq      ymm3, ymm3, ymm7
    
    vpmuldq     ymm4, ymm8, ymm12
    vpmuldq     ymm6, ymm10, ymm12
    vpmuldq     ymm5, ymm9, ymm12
    vpmuldq     ymm7, ymm11, ymm12
    
    vpmuldq     ymm8, ymm4, ymm14 //mul = reduce(a * sub)
    vpmuldq     ymm9, ymm5, ymm14                    
    vpmuldq     ymm10, ymm6, ymm14 
    vpmuldq     ymm11, ymm7, ymm14
    vpmuludq    ymm8, ymm8, ymm15  
    vpmuludq    ymm9, ymm9, ymm15 
    vpmuludq    ymm10, ymm10, ymm15                   
    vpmuludq    ymm11, ymm11, ymm15
    vpaddq      ymm4, ymm4, ymm8                      
    vpaddq      ymm5, ymm5, ymm9                     
    vpaddq      ymm6, ymm6, ymm10                                               
    vpaddq      ymm7, ymm7, ymm11  
    vpsrlq      ymm4, ymm4, 32                        
    vpsrlq      ymm5, ymm5, 32                       
    vpsrlq      ymm6, ymm6, 32                                      
    vpsrlq      ymm7, ymm7, 32 

    vmovdqa      [reg_p1+rax],      ymm0
    vmovdqa      [reg_p1+rax+8*4],  ymm1
    vmovdqa      [reg_p1+rax+8*8],  ymm2
    vmovdqa      [reg_p1+rax+8*12], ymm3
    vmovdqa      [reg_p1+rax+8*32], ymm4
    vmovdqa      [reg_p1+rax+8*36], ymm5
    vmovdqa      [reg_p1+rax+8*40], ymm6
    vmovdqa      [reg_p1+rax+8*44], ymm7

    add          rax, 8*16
    dec          r10
    jnz          loop_idgt5
    add          rax, 8*32
    inc          r8
    cmp          r8, 16
    jne          iround5loop
    
/* Round 6  */
    xor          rax, rax
    xor          r8, r8
iround6loop:
    mov          r10, 4
    vbroadcastss ymm12, DWORD PTR [reg_p3+4*r8]
loop_idgt6:
    vmovdqa      ymm0, [reg_p1+rax]
    vmovdqa      ymm1, [reg_p1+rax+8*64]
    vmovdqa      ymm2, [reg_p1+rax+8*4]
    vmovdqa      ymm3, [reg_p1+rax+8*68]
    vmovdqa      ymm4, [reg_p1+rax+8*8]
    vmovdqa      ymm5, [reg_p1+rax+8*72]
    vmovdqa      ymm6, [reg_p1+rax+8*12]
    vmovdqa      ymm7, [reg_p1+rax+8*76]
    
    vpsubd      ymm8, ymm0, ymm1
    vpsubd      ymm9, ymm2, ymm3
    vpsubd      ymm10, ymm4, ymm5
    vpsubd      ymm11, ymm6, ymm7
    
    vpaddd      ymm0, ymm0, ymm1
    vpaddd      ymm1, ymm2, ymm3
    vpaddd      ymm2, ymm4, ymm5
    vpaddd      ymm3, ymm6, ymm7

    vpmuldq     ymm4, ymm0, ymm13 // barr_reduce(add)
    vpmuldq     ymm5, ymm1, ymm13
    vpmuldq     ymm6, ymm2, ymm13
    vpmuldq     ymm7, ymm3, ymm13
    vpsrlq      ymm4, ymm4, PARAM_BARR_DIV
    vpsrlq      ymm5, ymm5, PARAM_BARR_DIV
    vpsrlq      ymm6, ymm6, PARAM_BARR_DIV
    vpsrlq      ymm7, ymm7, PARAM_BARR_DIV
    vpmuludq    ymm4, ymm4, ymm15
    vpmuludq    ymm5, ymm5, ymm15
    vpmuludq    ymm6, ymm6, ymm15
    vpmuludq    ymm7, ymm7, ymm15
    vpsubq      ymm0, ymm0, ymm4
    vpsubq      ymm1, ymm1, ymm5
    vpsubq      ymm2, ymm2, ymm6
    vpsubq      ymm3, ymm3, ymm7
    
    vpmuldq     ymm4, ymm8, ymm12
    vpmuldq     ymm6, ymm10, ymm12
    vpmuldq     ymm5, ymm9, ymm12
    vpmuldq     ymm7, ymm11, ymm12
    
    vpmuldq     ymm8, ymm4, ymm14 //mul = reduce(a * sub)
    vpmuldq     ymm9, ymm5, ymm14                    
    vpmuldq     ymm10, ymm6, ymm14 
    vpmuldq     ymm11, ymm7, ymm14
    vpmuludq    ymm8, ymm8, ymm15  
    vpmuludq    ymm9, ymm9, ymm15 
    vpmuludq    ymm10, ymm10, ymm15                   
    vpmuludq    ymm11, ymm11, ymm15
    vpaddq      ymm4, ymm4, ymm8                      
    vpaddq      ymm5, ymm5, ymm9                     
    vpaddq      ymm6, ymm6, ymm10                                               
    vpaddq      ymm7, ymm7, ymm11  
    vpsrlq      ymm4, ymm4, 32                        
    vpsrlq      ymm5, ymm5, 32                       
    vpsrlq      ymm6, ymm6, 32                                      
    vpsrlq      ymm7, ymm7, 32 

    vmovdqa     [reg_p1+rax],      ymm0
    vmovdqa     [reg_p1+rax+8*4],  ymm1
    vmovdqa     [reg_p1+rax+8*8],  ymm2
    vmovdqa     [reg_p1+rax+8*12], ymm3
    vmovdqa     [reg_p1+rax+8*64], ymm4
    vmovdqa     [reg_p1+rax+8*68], ymm5
    vmovdqa     [reg_p1+rax+8*72], ymm6
    vmovdqa     [reg_p1+rax+8*76], ymm7

    add         rax, 8*16
    dec         r10
    jnz         loop_idgt6
    add         rax, 8*64
    inc         r8
    cmp         r8, 8
    jne         iround6loop
    
    xor         rax, rax 
    vbroadcastss ymm12, DWORD PTR [reg_p3]
    mov         r10, 32
loop_idgt789:
/* Round 7  */
    vmovdqa      ymm0, [reg_p1+2*rax]
    vmovdqa      ymm1, [reg_p1+2*rax+8*128]
    vmovdqa      ymm2, [reg_p1+2*rax+8*256]
    vmovdqa      ymm3, [reg_p1+2*rax+8*384]
    vmovdqa      ymm4, [reg_p1+2*rax+8*512]
    vmovdqa      ymm5, [reg_p1+2*rax+8*640]
    vmovdqa      ymm6, [reg_p1+2*rax+8*768]
    vmovdqa      ymm7, [reg_p1+2*rax+8*896]
    
    vpsubd      ymm8, ymm0, ymm1
    vpsubd      ymm9, ymm2, ymm3
    vpsubd      ymm10, ymm4, ymm5
    vpsubd      ymm11, ymm6, ymm7
    
    vpaddd      ymm0, ymm0, ymm1
    vpaddd      ymm1, ymm2, ymm3
    vpaddd      ymm2, ymm4, ymm5
    vpaddd      ymm3, ymm6, ymm7

    vpmuldq     ymm4, ymm0, ymm13 // barr_reduce(add)
    vpmuldq     ymm5, ymm1, ymm13
    vpmuldq     ymm6, ymm2, ymm13
    vpmuldq     ymm7, ymm3, ymm13
    vpsrlq      ymm4, ymm4, PARAM_BARR_DIV
    vpsrlq      ymm5, ymm5, PARAM_BARR_DIV
    vpsrlq      ymm6, ymm6, PARAM_BARR_DIV
    vpsrlq      ymm7, ymm7, PARAM_BARR_DIV
    vpmuludq    ymm4, ymm4, ymm15
    vpmuludq    ymm5, ymm5, ymm15
    vpmuludq    ymm6, ymm6, ymm15
    vpmuludq    ymm7, ymm7, ymm15
    vpsubq      ymm0, ymm0, ymm4
    vpsubq      ymm1, ymm1, ymm5
    vpsubq      ymm2, ymm2, ymm6
    vpsubq      ymm3, ymm3, ymm7
    
    vpmuldq     ymm4, ymm8, ymm12
    vbroadcastss ymm8, DWORD PTR [reg_p3+4]
    vpmuldq     ymm5, ymm9, ymm8
    vbroadcastss ymm9, DWORD PTR [reg_p3+8]
    vpmuldq     ymm6, ymm10, ymm9
    vbroadcastss ymm10, DWORD PTR [reg_p3+12]
    vpmuldq     ymm7, ymm11, ymm10
    
    vpmuldq     ymm8, ymm4, ymm14 //mul = reduce(a * sub)
    vpmuldq     ymm9, ymm5, ymm14                    
    vpmuldq     ymm10, ymm6, ymm14 
    vpmuldq     ymm11, ymm7, ymm14
    vpmuludq    ymm8, ymm8, ymm15  
    vpmuludq    ymm9, ymm9, ymm15 
    vpmuludq    ymm10, ymm10, ymm15                   
    vpmuludq    ymm11, ymm11, ymm15
    vpaddq      ymm4, ymm4, ymm8                      
    vpaddq      ymm5, ymm5, ymm9                     
    vpaddq      ymm6, ymm6, ymm10                                               
    vpaddq      ymm7, ymm7, ymm11  
    vpsrlq      ymm4, ymm4, 32                        
    vpsrlq      ymm5, ymm5, 32                       
    vpsrlq      ymm6, ymm6, 32                                      
    vpsrlq      ymm7, ymm7, 32
    
    vpsubd      ymm8, ymm0, ymm1
    vpsubd      ymm9, ymm2, ymm3
    vpsubd      ymm10, ymm4, ymm5
    vpsubd      ymm11, ymm6, ymm7
    
    vpaddd      ymm0, ymm0, ymm1
    vpaddd      ymm1, ymm2, ymm3
    vpaddd      ymm2, ymm4, ymm5
    vpaddd      ymm3, ymm6, ymm7

    vpmuldq     ymm4, ymm0, ymm13 // barr_reduce(add)
    vpmuldq     ymm5, ymm1, ymm13
    vpmuldq     ymm6, ymm2, ymm13
    vpmuldq     ymm7, ymm3, ymm13
    vpsrlq      ymm4, ymm4, PARAM_BARR_DIV
    vpsrlq      ymm5, ymm5, PARAM_BARR_DIV
    vpsrlq      ymm6, ymm6, PARAM_BARR_DIV
    vpsrlq      ymm7, ymm7, PARAM_BARR_DIV
    vpmuludq    ymm4, ymm4, ymm15
    vpmuludq    ymm5, ymm5, ymm15
    vpmuludq    ymm6, ymm6, ymm15
    vpmuludq    ymm7, ymm7, ymm15
    vpsubq      ymm0, ymm0, ymm4
    vpsubq      ymm1, ymm1, ymm5
    vpsubq      ymm2, ymm2, ymm6
    vpsubq      ymm3, ymm3, ymm7
    
    vpmuldq     ymm4, ymm8, ymm12
    vpmuldq     ymm6, ymm10, ymm12
    vbroadcastss ymm8, DWORD PTR [reg_p3+4]
    vpmuldq     ymm5, ymm9, ymm8
    vpmuldq     ymm7, ymm11, ymm8
    
    vpmuldq     ymm8, ymm4, ymm14 //mul = reduce(a * sub)
    vpmuldq     ymm9, ymm5, ymm14                    
    vpmuldq     ymm10, ymm6, ymm14 
    vpmuldq     ymm11, ymm7, ymm14
    vpmuludq    ymm8, ymm8, ymm15  
    vpmuludq    ymm9, ymm9, ymm15 
    vpmuludq    ymm10, ymm10, ymm15                   
    vpmuludq    ymm11, ymm11, ymm15
    vpaddq      ymm4, ymm4, ymm8                      
    vpaddq      ymm5, ymm5, ymm9                     
    vpaddq      ymm6, ymm6, ymm10                                               
    vpaddq      ymm7, ymm7, ymm11  
    vpsrlq      ymm4, ymm4, 32                        
    vpsrlq      ymm5, ymm5, 32                       
    vpsrlq      ymm6, ymm6, 32                                      
    vpsrlq      ymm7, ymm7, 32
    
    vpsubd      ymm8, ymm0, ymm1
    vpsubd      ymm9, ymm2, ymm3
    vpsubd      ymm10, ymm4, ymm5
    vpsubd      ymm11, ymm6, ymm7
    
    vpaddd      ymm0, ymm0, ymm1
    vpaddd      ymm1, ymm2, ymm3
    vpaddd      ymm2, ymm4, ymm5
    vpaddd      ymm3, ymm6, ymm7

    vpmuldq     ymm4, ymm0, ymm13 // barr_reduce(add)
    vpmuldq     ymm5, ymm1, ymm13
    vpmuldq     ymm6, ymm2, ymm13
    vpmuldq     ymm7, ymm3, ymm13
    vpsrlq      ymm4, ymm4, PARAM_BARR_DIV
    vpsrlq      ymm5, ymm5, PARAM_BARR_DIV
    vpsrlq      ymm6, ymm6, PARAM_BARR_DIV
    vpsrlq      ymm7, ymm7, PARAM_BARR_DIV
    vpmuludq    ymm4, ymm4, ymm15
    vpmuludq    ymm5, ymm5, ymm15
    vpmuludq    ymm6, ymm6, ymm15
    vpmuludq    ymm7, ymm7, ymm15
    vpsubq      ymm0, ymm0, ymm4
    vpsubq      ymm1, ymm1, ymm5
    vpsubq      ymm2, ymm2, ymm6
    vpsubq      ymm3, ymm3, ymm7
    
    vpmuldq     ymm4, ymm8, ymm12
    vpmuldq     ymm6, ymm10, ymm12
    vpmuldq     ymm5, ymm9, ymm12
    vpmuldq     ymm7, ymm11, ymm12
    
    vpmuldq     ymm8, ymm4, ymm14 //mul = reduce(a * sub)
    vpmuldq     ymm9, ymm5, ymm14                    
    vpmuldq     ymm10, ymm6, ymm14 
    vpmuldq     ymm11, ymm7, ymm14
    vpmuludq    ymm8, ymm8, ymm15  
    vpmuludq    ymm9, ymm9, ymm15 
    vpmuludq    ymm10, ymm10, ymm15                   
    vpmuludq    ymm11, ymm11, ymm15
    vpaddq      ymm4, ymm4, ymm8                      
    vpaddq      ymm5, ymm5, ymm9                     
    vpaddq      ymm6, ymm6, ymm10                                               
    vpaddq      ymm7, ymm7, ymm11  
    vpsrlq      ymm4, ymm4, 32                        
    vpsrlq      ymm5, ymm5, 32                       
    vpsrlq      ymm6, ymm6, 32                                      
    vpsrlq      ymm7, ymm7, 32
    
    vmovdqa     [reg_p1+2*rax],       ymm0
    vmovdqa     [reg_p1+2*rax+8*128], ymm1
    vmovdqa     [reg_p1+2*rax+8*256], ymm2
    vmovdqa     [reg_p1+2*rax+8*384], ymm3
    vmovdqa     [reg_p1+2*rax+8*512], ymm4
    vmovdqa     [reg_p1+2*rax+8*640], ymm5
    vmovdqa     [reg_p1+2*rax+8*768], ymm6
    vmovdqa     [reg_p1+2*rax+8*896], ymm7

    add          rax, 4*4
    dec          r10
    jnz          loop_idgt789
    
    ret
