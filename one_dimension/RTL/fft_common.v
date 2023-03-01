//
// 采用基础算法描述 fft
//

// 这个模块使用了三个阶段来实现 FFT 算法。在第一阶段中，输入序列被分成偶数项和奇数项，
// 并存储在对应的寄存器中。在第二阶段中，通过递归调用FFT模块来计算偶数项和奇数项的
// DFT，并计算它们的加权和。在第三阶段中，将临时计算结果存储到输出寄存器中，以便最
// 终输出。在每个阶段中，都使用 Verilog 语言的条件语句和循环语句来实现相应的计算。
module fft_common (
    input clk, // 时钟信号
    input reset, // 复位信号
    input signed [31:0] x_r [0:N-1], // 输入实部
    input signed [31:0] x_i [0:N-1], // 输入虚部
    output reg signed [31:0] X_r [0:N-1], // 输出实部
    output reg signed [31:0] X_i [0:N-1] // 输出虚部
);

parameter N = 8; // 输入序列的长度，必须为2的幂次方
parameter log2N = $clog2(N);

complex_mult #(.N(N)) u[N/2]; // FFT中使用的复数乘法器
reg [31:0] x_even_r [0:N/2-1]; // 存储输入序列偶数项的实部
reg [31:0] x_even_i [0:N/2-1]; // 存储输入序列偶数项的虚部
reg [31:0] x_odd_r [0:N/2-1]; // 存储输入序列奇数项的实部
reg [31:0] x_odd_i [0:N/2-1]; // 存储输入序列奇数项的虚部
reg [31:0] W_r [0:N/2-1]; // 存储旋转因子的实部
reg [31:0] W_i [0:N/2-1]; // 存储旋转因子的虚部
reg [31:0] X_r_temp [0:N-1]; // 存储临时计算结果的实部
reg [31:0] X_i_temp [0:N-1]; // 存储临时计算结果的虚部
reg [31:0] X_r_final [0:N-1]; // 存储最终计算结果的实部
reg [31:0] X_i_final [0:N-1]; // 存储最终计算结果的虚部
reg [4:0] stage; // FFT的阶段

always @(posedge clk, posedge reset) begin
    if (reset) begin
        stage <= 0;
        x_even_r <= '{N{1'b0}};
        x_even_i <= '{N{1'b0}};
        x_odd_r <= '{N{1'b0}};
        x_odd_i <= '{N{1'b0}};
        W_r <= '{N/2{1'b0}};
        W_i <= '{N/2{1'b0}};
        X_r_temp <= '{N{1'b0}};
        X_i_temp <= '{N{1'b0}};
        X_r_final <= '{N{1'b0}};
        X_i_final <= '{N{1'b0}};
    end else begin
        case (stage)
            0: begin // 第一阶段：将输入序列分解为偶数项和奇数项
                for (i = 0; i < N/2; i = i + 1) begin
                    x_even_r[i] <= x_r[i*2];
                    x_even_i[i] <= x_i[i*2];
                    x_odd_r[i] <= x_r[i*2+1];
                    x_odd_i[i] <= x_i[i*2+1];
                end
                stage <= stage + 1;
            end
            1: begin // 第二阶段：递归计算偶数项和奇数项的DFT
                if (stage < log2N) begin
                    for (i = 0; i < N/2; i = i + 1) begin
                        u[i].a <= W_r[i];
                        u[i].b <= W_i[i];
                        u[i].x_in_r <= x_odd_r[i];
                        u[i].x_in_i <= x_odd_i[i];
                        u[i].y_in_r <= x_even_r[i];
                        u[i].y_in_i <= x_even_i[i];
                        u[i].clk <= clk;
                        u[i].reset <= reset;
                        X_r_temp[i] <= u[i].z_out_r;
                        X_i_temp[i] <= u[i].z_out_i;
                        X_r_temp[i+N/2] <= X_r_temp[i] - 2**(32-N)*W_i[i]*x_odd_i[i];
                        X_i_temp[i+N/2] <= X_i_temp[i] + 2**(32-N)*W_i[i]*x_odd_r[i];
                    end
                    stage <= stage + 1;
                end else begin
                    for (i = 0; i < N/2; i = i + 1) begin
                        X_r_temp[i] <= x_even_r[i] + x_odd_r[i];
                        X_i_temp[i] <= x_even_i[i] + x_odd_i[i];
                        X_r_temp[i+N/2] <= x_even_r[i] - x_odd_r[i];
                        X_i_temp[i+N/2] <= x_even_i[i] - x_odd_i[i];
                    end
                    stage <= stage + 1;
                end
            end
            2: begin // 第三阶段：将临时计算结果存储到输出寄存器中
                for (i = 0; i < N; i = i + 1) begin
                    X_r_final[i] <= X_r_temp[i];
                    X_i_final[i] <= X_i_temp[i];
                end
                stage <= 0;
            end
        endcase
    end
end

endmodule
