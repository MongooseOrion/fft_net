//
// ���û����㷨���� fft
//

// ���ģ��ʹ���������׶���ʵ�� FFT �㷨���ڵ�һ�׶��У��������б��ֳ�ż����������
// ���洢�ڶ�Ӧ�ļĴ����С��ڵڶ��׶��У�ͨ���ݹ����FFTģ��������ż������������
// DFT�����������ǵļ�Ȩ�͡��ڵ����׶��У�����ʱ�������洢������Ĵ����У��Ա���
// ���������ÿ���׶��У���ʹ�� Verilog ���Ե���������ѭ�������ʵ����Ӧ�ļ��㡣
module fft_common (
    input clk, // ʱ���ź�
    input reset, // ��λ�ź�
    input signed [31:0] x_r [0:N-1], // ����ʵ��
    input signed [31:0] x_i [0:N-1], // �����鲿
    output reg signed [31:0] X_r [0:N-1], // ���ʵ��
    output reg signed [31:0] X_i [0:N-1] // ����鲿
);

parameter N = 8; // �������еĳ��ȣ�����Ϊ2���ݴη�
parameter log2N = $clog2(N);

complex_mult #(.N(N)) u[N/2]; // FFT��ʹ�õĸ����˷���
reg [31:0] x_even_r [0:N/2-1]; // �洢��������ż�����ʵ��
reg [31:0] x_even_i [0:N/2-1]; // �洢��������ż������鲿
reg [31:0] x_odd_r [0:N/2-1]; // �洢���������������ʵ��
reg [31:0] x_odd_i [0:N/2-1]; // �洢����������������鲿
reg [31:0] W_r [0:N/2-1]; // �洢��ת���ӵ�ʵ��
reg [31:0] W_i [0:N/2-1]; // �洢��ת���ӵ��鲿
reg [31:0] X_r_temp [0:N-1]; // �洢��ʱ��������ʵ��
reg [31:0] X_i_temp [0:N-1]; // �洢��ʱ���������鲿
reg [31:0] X_r_final [0:N-1]; // �洢���ռ�������ʵ��
reg [31:0] X_i_final [0:N-1]; // �洢���ռ��������鲿
reg [4:0] stage; // FFT�Ľ׶�

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
            0: begin // ��һ�׶Σ����������зֽ�Ϊż�����������
                for (i = 0; i < N/2; i = i + 1) begin
                    x_even_r[i] <= x_r[i*2];
                    x_even_i[i] <= x_i[i*2];
                    x_odd_r[i] <= x_r[i*2+1];
                    x_odd_i[i] <= x_i[i*2+1];
                end
                stage <= stage + 1;
            end
            1: begin // �ڶ��׶Σ��ݹ����ż������������DFT
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
            2: begin // �����׶Σ�����ʱ�������洢������Ĵ�����
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
