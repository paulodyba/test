import sys
import time
import threading
from multiprocessing import Manager
from PyQt5.QtCore import QTimer
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QPushButton, QLabel


class WorkerThread(threading.Thread):
    def __init__(self, dev_id, shared_value, stop_event):
        super().__init__()
        self.dev_id = dev_id  # 唯一标识
        self.shared_value = shared_value  # 共享内存中的值，用来判断激活哪个线程
        self.stop_event = stop_event  # 停止线程的事件
        self._send_data = False  # 是否开始发送数据
        self.activated = False  # 记录线程是否已经激活

    def run(self):
        while not self.stop_event.is_set():
            if self.shared_value.value == self.dev_id and not self.activated:
                # 一旦激活后就开始发送数据，并且设置激活标志
                self.activated = True
                print(f"设备 {self.dev_id} 已激活，开始发送数据...")

            # 如果已经激活，始终发送数据
            if self.activated:
                print(f"设备 {self.dev_id} 正在持续发送数据...")

            # 线程一直在运行
            # print(f"设备 {self.dev_id} 正在运行。。。")
            time.sleep(1)  # 每秒检查一次

    def stop(self):
        """ 停止线程运行 """
        self.stop_event.set()  # 设置事件标志为已停止
        print(f"设备 {self.dev_id} 已停止")


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("线程控制示例")
        self.setGeometry(100, 100, 400, 300)

        # 使用 Manager 创建共享内存
        self.manager = Manager()
        self.shared_value = self.manager.Value('i', 0)  # 初始化共享内存，值为 0，表示第一个线程被激活
        self.shareda = 0

        self.threads = []
        self.stop_events = []  # 每个线程的停止事件
        self.thread_count = 6  # 需要创建6个线程
        self.status_label = QLabel("线程状态：", self)

        # 创建按钮
        self.start_threads_btn = QPushButton("启动所有线程", self)
        self.change_value_btn = QPushButton("改变共享值", self)

        # 创建六个停止按钮
        self.stop_buttons = []
        for i in range(self.thread_count):
            stop_button = QPushButton(f"停止设备 {i + 1}", self)
            stop_button.clicked.connect(self.create_stop_function(i))
            self.stop_buttons.append(stop_button)

        # 设置布局
        layout = QVBoxLayout()
        layout.addWidget(self.status_label)
        layout.addWidget(self.start_threads_btn)
        layout.addWidget(self.change_value_btn)
        for stop_button in self.stop_buttons:
            layout.addWidget(stop_button)

        self.setLayout(layout)

        # 连接信号与槽（信号槽分开管理）
        self.connect_signals()

    def connect_signals(self):
        """ 连接信号和槽 """
        self.start_threads_btn.clicked.connect(self.start_threads)
        self.change_value_btn.clicked.connect(self.change_shared_value)

    def start_threads(self):
        """ 启动所有线程 """
        for dev_id in range(self.thread_count):
            stop_event = threading.Event()  # 每个线程都有一个停止事件
            self.stop_events.append(stop_event)
            thread = WorkerThread(dev_id, self.shared_value, stop_event)
            thread.start()
            self.threads.append(thread)
        self.update_status(f"已启动 {self.thread_count} 个线程")

    def create_stop_function(self, dev_id):
        """ 创建停止按钮的槽函数 """

        def stop_thread():
            # 调用线程的 stop 方法
            self.threads[dev_id].stop()
            self.update_status(f"设备 {dev_id + 1} 已停止")

        return stop_thread

    def monitor_shared_value(self):
        """ 定时监测共享内存中的值，并根据值控制线程激活 """
        print(f"当前共享值: {self.shared_value.value}")
        # 打印当前状态
        for dev_id in range(self.thread_count):
            if self.shared_value.value == dev_id:
                print(f"设备 {dev_id + 1} 被激活")
            else:
                print(f"设备 {dev_id + 1} 继续运行")

    def change_shared_value(self):
        """ 改变共享内存中的值 """
        # 循环改变共享值，使其在 0 到 5 之间
        if self.shared_value.value < 5:
            self.shared_value.value += 1
        else:
            self.shared_value.value = 0
        self.update_status(f"共享值已更新为 {self.shared_value.value}")

    def update_status(self, message):
        """ 更新状态标签 """
        self.status_label.setText(message)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
