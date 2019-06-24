/*************************************************************************************
* qTESLA: an efficient post-quantum signature scheme based on the R-LWE problem
*
* Abstract: CDT constants for the Gaussian sampler
**************************************************************************************/

#ifndef CDTSAMP
#define CDTSAMP

#include <stdint.h>
#include "params.h"


// Sigma = 22.93, 64-bit precision

#define CDT_ROWS 207
#define CDT_COLS 2

static const int32_t cdt_v[CDT_ROWS*CDT_COLS] = {
    0x00000000L, 0x00000000L, // 0
    0x023A1B3FL, 0x4A499901L, // 1
    0x06AD3C4CL, 0x0CA08592L, // 2
    0x0B1D1E95L, 0x401E5DB9L, // 3
    0x0F879D85L, 0x73D5BFB7L, // 4
    0x13EA9C5CL, 0x2939948AL, // 5
    0x18440933L, 0x7FE9008DL, // 6
    0x1C91DFF1L, 0x48F0AE83L, // 7
    0x20D22D0FL, 0x100BC806L, // 8
    0x25031040L, 0x60F31377L, // 9
    0x2922BEEBL, 0x50B180CFL, // 10
    0x2D2F866AL, 0x1E289169L, // 11
    0x3127CE19L, 0x102CF7B2L, // 12
    0x350A1928L, 0x118E580DL, // 13
    0x38D5082CL, 0x6A7E620AL, // 14
    0x3C875A73L, 0x599D6D36L, // 15
    0x401FEF0EL, 0x33E6A3E9L, // 16
    0x439DC59EL, 0x183BDACEL, // 17
    0x46FFFEDAL, 0x27E0518BL, // 18
    0x4A45DCD3L, 0x174E5549L, // 19
    0x4D6EC2F3L, 0x49172E12L, // 20
    0x507A35C1L, 0x7D9AA338L, // 21
    0x5367DA64L, 0x752F8E31L, // 22
    0x563775EDL, 0x2DC9F137L, // 23
    0x58E8EC6BL, 0x2865CAFCL, // 24
    0x5B7C3FD0L, 0x5CCC8CBEL, // 25
    0x5DF18EA7L, 0x3326C087L, // 26
    0x6049129FL, 0x01DAE6B6L, // 27
    0x62831EF8L, 0x2B524213L, // 28
    0x64A01ED3L, 0x0A5D1038L, // 29
    0x66A09363L, 0x6544ED52L, // 30
    0x68851217L, 0x1F7909FBL, // 31
    0x6A4E42A8L, 0x589BF09CL, // 32
    0x6BFCDD30L, 0x162DC445L, // 33
    0x6D91A82DL, 0x7BCBF55CL, // 34
    0x6F0D7697L, 0x75D3528FL, // 35
    0x707125EDL, 0x13F82E79L, // 36
    0x71BD9C54L, 0x260C26C7L, // 37
    0x72F3C6C7L, 0x7D9C0191L, // 38
    0x74149755L, 0x04472E63L, // 39
    0x7521036DL, 0x21A138EAL, // 40
    0x761A0251L, 0x35015867L, // 41
    0x77008B94L, 0x30C0BD22L, // 42
    0x77D595B9L, 0x2DE3507FL, // 43
    0x789A14EEL, 0x19C5DB94L, // 44
    0x794EF9E2L, 0x6BE2990AL, // 45
    0x79F530BEL, 0x20A7F127L, // 46
    0x7A8DA031L, 0x08443399L, // 47
    0x7B1928A5L, 0x4D9D53CFL, // 48
    0x7B98A38CL, 0x72C68357L, // 49
    0x7C0CE2C7L, 0x5D698B25L, // 50
    0x7C76B02AL, 0x6EF32779L, // 51
    0x7CD6CD1DL, 0x09F74C79L, // 52
    0x7D2DF24DL, 0x5037123AL, // 53
    0x7D7CCF81L, 0x52E6CC5DL, // 54
    0x7DC40B76L, 0x6127DAEAL, // 55
    0x7E0443D9L, 0x16F11331L, // 56
    0x7E3E0D4BL, 0x48A00B90L, // 57
    0x7E71F37EL, 0x64E0EF47L, // 58
    0x7EA07957L, 0x6735C829L, // 59
    0x7ECA1921L, 0x78D7B202L, // 60
    0x7EEF44CBL, 0x639ED1AEL, // 61
    0x7F10662DL, 0x02BA119FL, // 62
    0x7F2DDF53L, 0x66EE6A14L, // 63
    0x7F480AD7L, 0x6F81453BL, // 64
    0x7F5F3C32L, 0x2587B359L, // 65
    0x7F73C018L, 0x34C60C54L, // 66
    0x7F85DCD8L, 0x6B4FC49DL, // 67
    0x7F95D2B9L, 0x3769ED08L, // 68
    0x7FA3DC55L, 0x2996B8DEL, // 69
    0x7FB02EFAL, 0x0EEEE30FL, // 70
    0x7FBAFB03L, 0x45D73B72L, // 71
    0x7FC46C34L, 0x7C8C59F2L, // 72
    0x7FCCAA10L, 0x15CAA326L, // 73
    0x7FD3D828L, 0x7BEA4849L, // 74
    0x7FDA1675L, 0x3608E7C2L, // 75
    0x7FDF819AL, 0x1D3DFF35L, // 76
    0x7FE43333L, 0x1952FF5FL, // 77
    0x7FE84217L, 0x5506F15AL, // 78
    0x7FEBC29AL, 0x61880546L, // 79
    0x7FEEC6C7L, 0x4786A8A8L, // 80
    0x7FF15E99L, 0x0A1CB795L, // 81
    0x7FF3982EL, 0x24C17DCCL, // 82
    0x7FF57FFAL, 0x11B43169L, // 83
    0x7FF720EFL, 0x69B7A428L, // 84
    0x7FF884ABL, 0x30B995E4L, // 85
    0x7FF9B396L, 0x651D9C1EL, // 86
    0x7FFAB50BL, 0x68EE9B1AL, // 87
    0x7FFB8F72L, 0x5D4208A6L, // 88
    0x7FFC485EL, 0x08AD19C4L, // 89
    0x7FFCE4A3L, 0x61DC95CCL, // 90
    0x7FFD6873L, 0x573AAF25L, // 91
    0x7FFDD76BL, 0x6C207ED1L, // 92
    0x7FFE34AAL, 0x43673438L, // 93
    0x7FFE82DEL, 0x2E535443L, // 94
    0x7FFEC454L, 0x55D51370L, // 95
    0x7FFEFB06L, 0x12FD6DC5L, // 96
    0x7FFF28A2L, 0x0A588B08L, // 97
    0x7FFF4E98L, 0x1CA2A14FL, // 98
    0x7FFF6E21L, 0x3E0B4535L, // 99
    0x7FFF8847L, 0x43F95CC4L, // 100
    0x7FFF9DEBL, 0x38044301L, // 101
    0x7FFFAFCBL, 0x3DA0CF24L, // 102
    0x7FFFBE88L, 0x16D5DC7CL, // 103
    0x7FFFCAA8L, 0x532DED04L, // 104
    0x7FFFD49EL, 0x330C43AAL, // 105
    0x7FFFDCC8L, 0x488C8B03L, // 106
    0x7FFFE376L, 0x5E2582C2L, // 107
    0x7FFFE8EBL, 0x2A699905L, // 108
    0x7FFFED5DL, 0x5773C7A7L, // 109
    0x7FFFF0FBL, 0x63D3499FL, // 110
    0x7FFFF3EBL, 0x621D490AL, // 111
    0x7FFFF64DL, 0x1BAFE266L, // 112
    0x7FFFF83AL, 0x1AA50219L, // 113
    0x7FFFF9C8L, 0x1E74DD87L, // 114
    0x7FFFFB08L, 0x7E5630D3L, // 115
    0x7FFFFC0AL, 0x7C050D38L, // 116
    0x7FFFFCDAL, 0x093EEF3BL, // 117
    0x7FFFFD80L, 0x01F3172BL, // 118
    0x7FFFFE04L, 0x5CDFCE2EL, // 119
    0x7FFFFE6EL, 0x54177CDFL, // 120
    0x7FFFFEC3L, 0x06B266A3L, // 121
    0x7FFFFF06L, 0x14C2B342L, // 122
    0x7FFFFF3BL, 0x367771F9L, // 123
    0x7FFFFF65L, 0x4F37BDD3L, // 124
    0x7FFFFF86L, 0x7D6081B5L, // 125
    0x7FFFFFA1L, 0x2734F6F5L, // 126
    0x7FFFFFB6L, 0x057B565CL, // 127
    0x7FFFFFC6L, 0x2C2BD768L, // 128
    0x7FFFFFD3L, 0x118798A8L, // 129
    0x7FFFFFDDL, 0x13DF050CL, // 130
    0x7FFFFFE4L, 0x7E436700L, // 131
    0x7FFFFFEBL, 0x0C554F26L, // 132
    0x7FFFFFEFL, 0x6D58FEBAL, // 133
    0x7FFFFFF3L, 0x46B2EA4DL, // 134
    0x7FFFFFF6L, 0x35E875C6L, // 135
    0x7FFFFFF8L, 0x523C11B9L, // 136
    0x7FFFFFFAL, 0x2DF7BE14L, // 137
    0x7FFFFFFBL, 0x577585A6L, // 138
    0x7FFFFFFCL, 0x59F2AC82L, // 139
    0x7FFFFFFDL, 0x3E37F0C9L, // 140
    0x7FFFFFFEL, 0x0B1F4CF2L, // 141
    0x7FFFFFFEL, 0x45FE12ACL, // 142
    0x7FFFFFFEL, 0x72F8E740L, // 143
    0x7FFFFFFFL, 0x154618FFL, // 144
    0x7FFFFFFFL, 0x2F61E68CL, // 145
    0x7FFFFFFFL, 0x43379BB6L, // 146
    0x7FFFFFFFL, 0x5241D483L, // 147
    0x7FFFFFFFL, 0x5DA3C063L, // 148
    0x7FFFFFFFL, 0x663CDF59L, // 149
    0x7FFFFFFFL, 0x6CB865F1L, // 150
    0x7FFFFFFFL, 0x71993691L, // 151
    0x7FFFFFFFL, 0x75432D5CL, // 152
    0x7FFFFFFFL, 0x780253E4L, // 153
    0x7FFFFFFFL, 0x7A10727DL, // 154
    0x7FFFFFFFL, 0x7B995BC9L, // 155
    0x7FFFFFFFL, 0x7CBE3B28L, // 156
    0x7FFFFFFFL, 0x7D981EEFL, // 157
    0x7FFFFFFFL, 0x7E39EAD2L, // 158
    0x7FFFFFFFL, 0x7EB1D52CL, // 159
    0x7FFFFFFFL, 0x7F0A8A07L, // 160
    0x7FFFFFFFL, 0x7F4C08CCL, // 161
    0x7FFFFFFFL, 0x7F7C4CC9L, // 162
    0x7FFFFFFFL, 0x7F9FCD06L, // 163
    0x7FFFFFFFL, 0x7FB9DD06L, // 164
    0x7FFFFFFFL, 0x7FCCF5DEL, // 165
    0x7FFFFFFFL, 0x7FDAED50L, // 166
    0x7FFFFFFFL, 0x7FE51F3EL, // 167
    0x7FFFFFFFL, 0x7FEC8CC3L, // 168
    0x7FFFFFFFL, 0x7FF1F385L, // 169
    0x7FFFFFFFL, 0x7FF5DF23L, // 170
    0x7FFFFFFFL, 0x7FF8B62FL, // 171
    0x7FFFFFFFL, 0x7FFAC3DFL, // 172
    0x7FFFFFFFL, 0x7FFC3F40L, // 173
    0x7FFFFFFFL, 0x7FFD5084L, // 174
    0x7FFFFFFFL, 0x7FFE14FBL, // 175
    0x7FFFFFFFL, 0x7FFEA1F4L, // 176
    0x7FFFFFFFL, 0x7FFF06ECL, // 177
    0x7FFFFFFFL, 0x7FFF4F19L, // 178
    0x7FFFFFFFL, 0x7FFF8298L, // 179
    0x7FFFFFFFL, 0x7FFFA744L, // 180
    0x7FFFFFFFL, 0x7FFFC155L, // 181
    0x7FFFFFFFL, 0x7FFFD3D3L, // 182
    0x7FFFFFFFL, 0x7FFFE0EBL, // 183
    0x7FFFFFFFL, 0x7FFFEA2CL, // 184
    0x7FFFFFFFL, 0x7FFFF0B3L, // 185
    0x7FFFFFFFL, 0x7FFFF54CL, // 186
    0x7FFFFFFFL, 0x7FFFF886L, // 187
    0x7FFFFFFFL, 0x7FFFFACAL, // 188
    0x7FFFFFFFL, 0x7FFFFC60L, // 189
    0x7FFFFFFFL, 0x7FFFFD7CL, // 190
    0x7FFFFFFFL, 0x7FFFFE42L, // 191
    0x7FFFFFFFL, 0x7FFFFECBL, // 192
    0x7FFFFFFFL, 0x7FFFFF2BL, // 193
    0x7FFFFFFFL, 0x7FFFFF6DL, // 194
    0x7FFFFFFFL, 0x7FFFFF9BL, // 195
    0x7FFFFFFFL, 0x7FFFFFBBL, // 196
    0x7FFFFFFFL, 0x7FFFFFD1L, // 197
    0x7FFFFFFFL, 0x7FFFFFE0L, // 198
    0x7FFFFFFFL, 0x7FFFFFEAL, // 199
    0x7FFFFFFFL, 0x7FFFFFF1L, // 200
    0x7FFFFFFFL, 0x7FFFFFF6L, // 201
    0x7FFFFFFFL, 0x7FFFFFF9L, // 202
    0x7FFFFFFFL, 0x7FFFFFFCL, // 203
    0x7FFFFFFFL, 0x7FFFFFFDL, // 204
    0x7FFFFFFFL, 0x7FFFFFFEL, // 205
    0x7FFFFFFFL, 0x7FFFFFFFL, // 206
}; // cdt_v

// memory requirements:
//      512 samples:  8628 bytes
//      256 samples:  5556 bytes
//      128 samples:  4020 bytes
//       64 samples:  3252 bytes
//       32 samples:  2868 bytes
// table alone: 1656 bytes

#endif 