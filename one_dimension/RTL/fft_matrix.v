//////////////////////////////
// һά����Ҷ�任ģ�����
//////////////////////////////

// ��ģ������8λʮ���Ƶ����ݣ�����ʵ�����鲿����������
module fft_matrix(
  input clk,
  input rst,
  input [255:0] x_r, // ����ʵ��
  input [255:0] x_i, // �����鲿
  output reg [255:0] y_r, // ���ʵ��
  output reg [255:0] y_i // ����鲿
);

  // ���峣��
  parameter N = 256; // ���ݵ���
  parameter LOG2N = 8; // N�Ķ���
  parameter PI = 3.141592653589793;

  // �������
  reg signed [7:0] X [0:N-1][0:1]; // �洢��������
  reg signed [7:0] Y [0:N-1][0:1]; // �洢�м���
  reg signed [7:0] W_r [0:N/2-1], W_i [0:N/2-1]; // �洢��ת����
  reg signed [7:0] T [0:N-1][0:N-1]; // �洢����T

  // ��ʼ����ת����
  initial begin
    for (int k = 0; k < N/2; k = k + 1) begin
      W_r[k] = $signed(8'hFF & (cos(2 * PI * k / N) * 128));
      W_i[k] = $signed(8'hFF & (-sin(2 * PI * k / N) * 128));
    end
  end

  // ��ʼ����������
  always @* begin
    for (int k = 0; k < N; k = k + 1) begin
      X[k][0] = x_r[k];
      X[k][1] = x_i[k];
    end
  end

  // ���ɾ���T
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

  // FFTģ��
  always @* begin
    Y = X * T;

    for (int k = 0; k < N; k = k + 1) begin
      y_r[k] = Y[k][0];
      y_i[k] = Y[k][1];
    end
  end

endmodule
