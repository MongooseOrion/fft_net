//////////////////////////////
// 一维傅里叶变换模块代码
//////////////////////////////

// 该模块允许8位十进制的数据（包含实部和虚部）输入运算
module fft_matrix(
  input clk,
  input rst,
  input [255:0] x_r, // 输入实部
  input [255:0] x_i, // 输入虚部
  output reg [255:0] y_r, // 输出实部
  output reg [255:0] y_i // 输出虚部
);

  // 定义常数
  parameter N = 256; // 数据点数
  parameter LOG2N = 8; // N的对数
  parameter PI = 3.141592653589793;

  // 定义矩阵
  reg signed [7:0] X [0:N-1][0:1]; // 存储输入数据
  reg signed [7:0] Y [0:N-1][0:1]; // 存储中间结果
  reg signed [7:0] W_r [0:N/2-1], W_i [0:N/2-1]; // 存储旋转因子
  reg signed [7:0] T [0:N-1][0:N-1]; // 存储矩阵T

  // 初始化旋转因子
  initial begin
    for (int k = 0; k < N/2; k = k + 1) begin
      W_r[k] = $signed(8'hFF & (cos(2 * PI * k / N) * 128));
      W_i[k] = $signed(8'hFF & (-sin(2 * PI * k / N) * 128));
    end
  end

  // 初始化输入数据
  always @* begin
    for (int k = 0; k < N; k = k + 1) begin
      X[k][0] = x_r[k];
      X[k][1] = x_i[k];
    end
  end

  // 生成矩阵T
  generate
    for (int i = 0; i < LOG2N; i = i + 1) begin
      assign T[0][0] = 1;
      for (int j = 1; j < N; j = j + 1) begin
        assign T[0][j] = 0;
      end

      for (int j = 1; j < N; j = j + 1) begin
        reg signed [7:0] arg;
        arg = $signed(8'hFF & (j >> i));
        arg = arg & (N-1);
        T[j][0] = $signed(8'hFF & cos(2 * PI * arg / N));
        T[j][1] = $signed(8'hFF & (-sin(2 * PI * arg / N)));
      end
    end
  endgenerate

  // FFT模块
  always @* begin
    Y = X * T;

    for (int k = 0; k < N; k = k + 1) begin
      y_r[k] = Y[k][0];
      y_i[k] = Y[k][1];
    end
  end

endmodule
